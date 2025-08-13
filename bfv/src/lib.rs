//! Implementation of BFV https://eprint.iacr.org/2012/144.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::ops;

use arith::{Ring, RingParam, Rq, R};

// error deviation for the Gaussian(Normal) distribution
// sigma=3.2 from: https://eprint.iacr.org/2022/162.pdf page 5
const ERR_SIGMA: f64 = 3.2;

#[derive(Clone, Copy, Debug)]
pub struct Param {
    ring: RingParam,
    t: u64,
    p: u64,
}
impl Param {
    // returns the plaintext param
    pub fn pt(&self) -> RingParam {
        RingParam {
            q: self.t,
            n: self.ring.n,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SecretKey(Rq);

#[derive(Clone, Debug)]
pub struct PublicKey(Rq, Rq);

/// Relinearization key
#[derive(Clone, Debug)]
pub struct RLK(Rq, Rq);

// RLWE ciphertext
#[derive(Clone, Debug)]
pub struct RLWE(Rq, Rq);

impl RLWE {
    fn add(lhs: Self, rhs: Self) -> Self {
        RLWE(lhs.0 + rhs.0, lhs.1 + rhs.1)
    }
    pub fn remodule(&self, p: u64) -> RLWE {
        let x = self.0.remodule(p);
        let y = self.1.remodule(p);
        RLWE(x, y)
    }

    fn tensor(t: u64, a: &Self, b: &Self) -> (Rq, Rq, Rq) {
        let (q, n) = (a.0.param.q, a.0.param.n);
        // expand Q->PQ // TODO rm

        // get the coefficients in Z, ie. interpret a,b \in R (instead of R_q)
        let a0: R = a.0.clone().to_r(); // TODO rm clone()
        let a1: R = a.1.clone().to_r();
        let b0: R = b.0.clone().to_r();
        let b1: R = b.1.clone().to_r();

        // tensor (\in R) (2021-204 p.9)
        // NOTE: here can use *, but at first versions want to make it explicit
        // that we're using the naive mul. TODO use *.
        use arith::ring_n::naive_mul;
        let c0: Vec<i64> = naive_mul(&a0, &b0);
        let c1_l: Vec<i64> = naive_mul(&a0, &b1);
        let c1_r = naive_mul(&a1, &b0);
        let c1: Vec<i64> = itertools::zip_eq(c1_l, c1_r).map(|(l, r)| l + r).collect();
        let c2: Vec<i64> = naive_mul(&a1, &b1);

        // scale down, then reduce module Q, so result is \in R_q
        let c0: Rq = arith::ring_n::mul_div_round(q, n, c0, t, q);
        let c1: Rq = arith::ring_n::mul_div_round(q, n, c1, t, q);
        let c2: Rq = arith::ring_n::mul_div_round(q, n, c2, t, q);

        (c0, c1, c2)
    }
    /// ciphertext multiplication
    fn mul(t: u64, rlk: &RLK, a: &Self, b: &Self) -> Self {
        let (c0, c1, c2) = Self::tensor(t, a, b);
        BFV::relinearize_204(&rlk, &c0, &c1, &c2)
    }
}
// naive mul in the ring Rq, reusing the ring_n::naive_mul and then applying mod(X^N +1)
fn tmp_naive_mul(a: Rq, b: Rq) -> Rq {
    Rq::from_vec_i64(
        &a.param.clone(),
        arith::ring_n::naive_mul(&a.to_r(), &b.to_r()),
    )
}

impl ops::Add<RLWE> for RLWE {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::add(self, rhs)
    }
}

impl ops::Add<&Rq> for &RLWE {
    type Output = RLWE;
    fn add(self, rhs: &Rq) -> Self::Output {
        BFV::add_const(self, rhs)
    }
}

pub struct BFV {}

impl BFV {
    // const DELTA: u64 = Q / T; // floor

    /// generate a new key pair (privK, pubK)
    pub fn new_key(mut rng: impl Rng, param: &Param) -> Result<(SecretKey, PublicKey)> {
        // WIP: review probabilities

        // let Xi_key = Uniform::new(-1_f64, 1_f64);
        let Xi_key = Uniform::new(0_u64, 2_u64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        // secret key
        // let mut s = Rq::rand_f64(&mut rng, Xi_key)?;
        let mut s = Rq::rand_u64(&mut rng, Xi_key, &param.ring)?;
        // since s is going to be multiplied by other Rq elements, already
        // compute its NTT
        s.compute_evals();

        // pk = (-a * s + e, a)
        let a = Rq::rand_u64(&mut rng, Uniform::new(0_u64, param.ring.q), &param.ring)?;
        let e = Rq::rand_f64(&mut rng, Xi_err, &param.ring)?;
        let pk: PublicKey = PublicKey(&(&(-a.clone()) * &s) + &e, a.clone()); // TODO rm clones
        Ok((SecretKey(s), pk))
    }

    // note: m is modulus t
    pub fn encrypt(mut rng: impl Rng, param: &Param, pk: &PublicKey, m: &Rq) -> Result<RLWE> {
        // assert param & inputs
        debug_assert_eq!(param.ring, pk.0.param);
        debug_assert_eq!(param.t, m.param.q);
        debug_assert_eq!(param.ring.n, m.param.n);

        let Xi_key = Uniform::new(-1_f64, 1_f64);
        // let Xi_key = Uniform::new(0_u64, 2_u64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u = Rq::rand_f64(&mut rng, Xi_key, &param.ring)?;
        // let u = Rq::rand_u64(&mut rng, Xi_key)?;
        let e_1 = Rq::rand_f64(&mut rng, Xi_err, &param.ring)?;
        let e_2 = Rq::rand_f64(&mut rng, Xi_err, &param.ring)?;

        // migrate m's coeffs to the bigger modulus Q (from T)
        let m = m.remodule(param.ring.q);
        let c0 = &pk.0 * &u + e_1 + m * (param.ring.q / param.t); // floor(q/t)=DELTA
        let c1 = &pk.1 * &u + e_2;
        Ok(RLWE(c0, c1))
    }

    pub fn decrypt(param: &Param, sk: &SecretKey, c: &RLWE) -> Rq {
        debug_assert_eq!(param.ring, sk.0.param);
        debug_assert_eq!(param.ring.q, c.0.param.q);
        debug_assert_eq!(param.ring.n, c.0.param.n);

        let cs: Rq = &c.0 + &(&c.1 * &sk.0); // done in mod q

        // same but with naive_mul:
        // let c1s = arith::ring_n::naive_mul(&c.1.to_r(), &sk.0.to_r());
        // let c1s = Rq::from_vec_i64(c1s);
        // let cs = c.0 + c1s;

        let r: Rq = cs.mul_div_round(param.t, param.ring.q);
        r.remodule(param.t)
    }

    fn add_const(c: &RLWE, m: &Rq) -> RLWE {
        let q = c.0.param.q;
        let t = m.param.q;

        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule(c.0.param.q);
        // TODO rm clones
        RLWE(c.0.clone() + m * (q / t), c.1.clone()) // floor(q/t)=DELTA
    }
    fn mul_const(rlk: &RLK, c: &RLWE, m: &Rq) -> RLWE {
        // let pq = rlk.0.q;
        let q = c.0.param.q;
        let t = m.param.q;

        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule(q);

        // encrypt m*Delta without noise, and then perform normal ciphertext multiplication
        let md = RLWE(m * (q / t), Rq::zero(&c.0.param)); // floor(q/t)=DELTA
        RLWE::mul(t, &rlk, &c, &md)
    }

    fn rlk_key(mut rng: impl Rng, param: &Param, s: &SecretKey) -> Result<RLK> {
        let pq = param.p * param.ring.q;
        let rlk_param = RingParam {
            q: pq,
            n: param.ring.n,
        };

        // TODO review using Xi' instead of Xi
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;
        // let Xi_err = Normal::new(0_f64, 0.0)?;
        let s = s.0.remodule(pq);
        let a = Rq::rand_u64(&mut rng, Uniform::new(0_u64, pq), &rlk_param)?;
        let e = Rq::rand_f64(&mut rng, Xi_err, &rlk_param)?;

        // let rlk: RLK<PQ, N> = RLK::<PQ, N>(-(&a * &s + e) + (s * s) * P, a.clone());
        // TODO rm clones
        let rlk: RLK = RLK(
            -(tmp_naive_mul(a.clone(), s.clone()) + e)
                + tmp_naive_mul(s.clone(), s.clone()) * param.p,
            a.clone(),
        );

        Ok(rlk)
    }

    fn relinearize(rlk: &RLK, c0: &Rq, c1: &Rq, c2: &Rq) -> RLWE {
        let pq = rlk.0.param.q;
        let param = c0.param;
        let q = param.q;
        let p = pq / q;

        let c2rlk0: Vec<f64> = (c2.clone().to_r() * rlk.0.clone().to_r())
            .coeffs()
            .iter()
            .map(|e| (*e as f64 / p as f64).round())
            .collect();

        let c2rlk1: Vec<f64> = (c2.clone().to_r() * rlk.1.clone().to_r()) // TODO rm clones
            .coeffs()
            .iter()
            .map(|e| (*e as f64 / p as f64).round())
            .collect();

        let r0 = Rq::from_vec_f64(&param, c2rlk0);
        let r1 = Rq::from_vec_f64(&param, c2rlk1);

        let res = RLWE(c0 + &r0, c1 + &r1);
        res
    }
    fn relinearize_204(rlk: &RLK, c0: &Rq, c1: &Rq, c2: &Rq) -> RLWE {
        let pq = rlk.0.param.q;
        let q = c0.param.q;
        let p = pq / q;
        let n = c0.param.n;
        // TODO (in debug) check that all Ns match

        // let c2rlk0: Rq<PQ, N> = c2.remodule::<PQ>() * rlk.0.remodule::<PQ>();
        // let c2rlk1: Rq<PQ, N> = c2.remodule::<PQ>() * rlk.1.remodule::<PQ>();
        // let r0: Rq<Q, N> = c2rlk0.mul_div_round(1, P).remodule::<Q>();
        // let r1: Rq<Q, N> = c2rlk1.mul_div_round(1, P).remodule::<Q>();

        use arith::ring_n::naive_mul;
        let c2rlk0: Vec<i64> = naive_mul(&c2.clone().to_r(), &rlk.0.clone().to_r()); // TODO rm clones
        let c2rlk1: Vec<i64> = naive_mul(&c2.clone().to_r(), &rlk.1.clone().to_r());
        let r0: Rq = arith::ring_n::mul_div_round(q, n, c2rlk0, 1, p);
        let r1: Rq = arith::ring_n::mul_div_round(q, n, c2rlk1, 1, p);

        let res = RLWE(c0 + &r0, c1 + &r1);
        res
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1, // q prime, and 2^q + 1 shape
                n: 512,
            },
            t: 32, // plaintext modulus
            p: 0,  // unused in this test
        };

        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let (sk, pk) = BFV::new_key(&mut rng, &param)?;

            let msg_dist = Uniform::new(0_u64, param.t);
            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

            let c = BFV::encrypt(&mut rng, &param, &pk, &m)?;
            let m_recovered = BFV::decrypt(&param, &sk, &c);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1, // q prime, and 2^q + 1 shape
                n: 128,
            },
            t: 32, // plaintext modulus
            p: 0,  // unused in this test
        };

        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let (sk, pk) = BFV::new_key(&mut rng, &param)?;

            let msg_dist = Uniform::new(0_u64, param.t);
            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

            let c1 = BFV::encrypt(&mut rng, &param, &pk, &m1)?;
            let c2 = BFV::encrypt(&mut rng, &param, &pk, &m2)?;

            let c3 = c1 + c2;

            let m3_recovered = BFV::decrypt(&param, &sk, &c3);

            assert_eq!(m1 + m2, m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_constant_add_mul() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1; // q prime, and 2^q + 1 shape
        let param = Param {
            ring: RingParam { q, n: 16 },
            t: 8, // plaintext modulus
            p: q * q,
        };

        let mut rng = rand::thread_rng();

        let (sk, pk) = BFV::new_key(&mut rng, &param)?;

        let msg_dist = Uniform::new(0_u64, param.t);
        let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
        let m2_const = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
        let c1 = BFV::encrypt(&mut rng, &param, &pk, &m1)?;

        let c3_add = &c1 + &m2_const;

        let m3_add_recovered = BFV::decrypt(&param, &sk, &c3_add);
        assert_eq!(&m1 + &m2_const, m3_add_recovered);

        // test multiplication of a ciphertext by a constant
        let rlk = BFV::rlk_key(&mut rng, &param, &sk)?;

        let c3_mul = BFV::mul_const(&rlk, &c1, &m2_const);

        let m3_mul_recovered = BFV::decrypt(&param, &sk, &c3_mul);
        assert_eq!(
            (m1.to_r() * m2_const.to_r()).to_rq(param.t).coeffs(),
            m3_mul_recovered.coeffs()
        );

        Ok(())
    }

    /*
        // TMP WIP
        #[test]
        #[ignore]
        fn test_param() -> Result<()> {
            const Q: u64 = 2u64.pow(16) + 1; // q prime, and 2^q + 1 shape
            const N: usize = 32;
            const T: u64 = 8; // plaintext modulus

            const P: u64 = Q * Q;
            const PQ: u64 = P * Q;
            const DELTA: u64 = Q / T; // floor

            let mut rng = rand::thread_rng();

            let Xi_key = Uniform::new(0_f64, 1_f64);
            let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

            let s = Rq::rand_f64(&mut rng, Xi_key)?;
            let e = Rq::rand_f64(&mut rng, Xi_err)?;
            let u = Rq::rand_f64(&mut rng, Xi_key)?;
            let e_0 = Rq::rand_f64(&mut rng, Xi_err)?;
            let e_1 = Rq::rand_f64(&mut rng, Xi_err)?;
            let m = Rq::rand_u64(&mut rng, Uniform::new(0_u64, T))?;

            // v_fresh
            let v: Rq<Q, N> = u * e + e_1 * s + e_0;

            let q: f64 = Q as f64;
            let t: f64 = T as f64;
            let n: f64 = N as f64;
            let delta: f64 = DELTA as f64;

            // r_t(q)/t should be equal to q/t-Δ
            assert_eq!(
                // r_t(q)/t, where r_t(q)=q mod t
                (q % t) / t,
                // Δt/Q = q - r_t(Q)/Q, so r_t(Q)=q - Δt
                (q / t) - delta
            );
            let rt: f64 = (q % t) / t;
            dbg!(&rt);

            dbg!(v.infinity_norm());
            let bound: f64 = (q / (2_f64 * t)) - (rt / 2_f64);
            dbg!(bound);
            assert!((v.infinity_norm() as f64) < bound);
            let max_v_infnorm = bound - 1.0;

            // addition noise
            let v_add: Rq<Q, N> = v + v + u * rt;
            let v_add: Rq<Q, N> = v_add + v_add + u * rt;
            assert!((v_add.infinity_norm() as f64) < bound);

            // multiplication noise
            let (_, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;
            let c = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m.remodule::<T>())?;
            let b_key: f64 = 1_f64;
            // ef: expansion factor
            let ef: f64 = 2.0 * n.sqrt();
            let bound: f64 = ((ef * t) / 2.0)
                * ((2.0 * max_v_infnorm * max_v_infnorm) / q
                    + (4.0 + ef * b_key) * (max_v_infnorm + max_v_infnorm)
                    + rt * (ef * b_key + 5.0))
                + (1.0 + ef * b_key + ef * ef * b_key * b_key) / 2.0;
            dbg!(&bound);

            let k: Vec<f64> = (c.0 + c.1 * s - m * delta - v)
                .coeffs()
                .iter()
                .map(|e_i| e_i.0 as f64 / q)
                .collect();
            let k = Rq::from_vec_f64(k);
            let v_tensor_0 = (v * v)
                .coeffs()
                .iter()
                .map(|e_i| (e_i.0 as f64 * t) / q)
                .collect::<Vec<f64>>();
            let v_tensor_0 = Rq::from_vec_f64(v_tensor_0);
            let v_tensor_1 = ((m * v) + (m * v))
                .coeffs()
                .iter()
                .map(|e_i| (e_i.0 as f64 * t * delta) / q)
                .collect::<Vec<f64>>();
            let v_tensor_1 = Rq::from_vec_f64(v_tensor_1);
            let v_tensor_2: Rq<Q, N> = (v * k + v * k) * t;
            let rm: f64 = (ef * t) / 2.0;
            let rm: Rq<Q, N> = Rq::from_vec_f64(vec![rm; N]);
            let v_tensor_3: Rq<Q, N> = (m * k
                + m * k
                + rm
                + Rq::from_vec_f64(
                    ((m * m) * DELTA)
                        .coeffs()
                        .iter()
                        .map(|e_i| e_i.0 as f64 / q)
                        .collect::<Vec<f64>>(),
                ))
                * rt;
            let v_tensor = v_tensor_0 + v_tensor_1 + v_tensor_2 - v_tensor_3;

            let v_r = (1.0 + ef * b_key + ef * ef * b_key * b_key) / 2.0;
            let v_mult_norm = v_tensor.infinity_norm() as f64 + v_r;
            dbg!(&v_mult_norm);
            dbg!(&bound);
            assert!(v_mult_norm < bound);

            // let m1 = Rq::<T, N>::zero();
            // let m2 = Rq::<T, N>::zero();
            // let (_, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;
            // let c1 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m1)?;
            // let c2 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m2)?;
            // let (c_a, c_b, c_c) = RLWE::tensor::<PQ, T>(&c1, &c2);
            // dbg!(&c_a.infinity_norm());
            // dbg!(&c_b.infinity_norm());
            // dbg!(&c_c.infinity_norm());
            // assert!((c_a.infinity_norm() as f64) < bound);
            // assert!((c_b.infinity_norm() as f64) < bound);
            // assert!((c_c.infinity_norm() as f64) < bound);
            // WIP

            Ok(())
        }
    */

    #[test]
    fn test_tensor() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1; // q prime, and 2^q + 1 shape
        let param = Param {
            ring: RingParam { q, n: 16 },
            t: 2, // plaintext modulus
            p: q * q,
        };
        let mut rng = rand::thread_rng();

        let msg_dist = Uniform::new(0_u64, param.t);
        for _ in 0..1_000 {
            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

            test_tensor_opt(&mut rng, &param, m1, m2)?;
        }

        Ok(())
    }
    fn test_tensor_opt(mut rng: impl Rng, param: &Param, m1: Rq, m2: Rq) -> Result<()> {
        let (sk, pk) = BFV::new_key(&mut rng, &param)?;

        let c1 = BFV::encrypt(&mut rng, &param, &pk, &m1)?;
        let c2 = BFV::encrypt(&mut rng, &param, &pk, &m2)?;

        let (c_a, c_b, c_c) = RLWE::tensor(param.t, &c1, &c2);
        // let (c_a, c_b, c_c) = RLWE::tensor_new::<PQ, T>(&c1, &c2);

        // decrypt non-relinearized mul result
        let m3: Rq = c_a + &c_b * &sk.0 + &c_c * &(&sk.0 * &sk.0);

        // let m3: Rq<Q, N> = c_a
        //     + Rq::from_vec_i64(arith::ring_n::naive_mul(&c_b.to_r(), &sk.0.to_r()))
        //     + Rq::from_vec_i64(arith::ring_n::naive_mul(
        //         &c_c.to_r(),
        //         &R::<N>::from_vec(arith::ring_n::naive_mul(&sk.0.to_r(), &sk.0.to_r())),
        //     ));
        let m3: Rq = m3.mul_div_round(param.t, param.ring.q); // descale
        let m3 = m3.remodule(param.t);

        let naive = (m1.clone().to_r() * m2.clone().to_r()).to_rq(param.t); // TODO rm clones
        assert_eq!(
            m3.coeffs().to_vec(),
            naive.coeffs().to_vec(),
            "\n\nfor testing:\nlet m1 = Rq::<T, N>::from_vec_u64(vec!{:?});\nlet m2 = Rq::<T, N>::from_vec_u64(vec!{:?});\n",
            m1.coeffs(),
            m2.coeffs()
        );

        Ok(())
    }

    #[test]
    fn test_mul_relin() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1; // q prime, and 2^q + 1 shape
        let param = Param {
            ring: RingParam { q, n: 16 },
            t: 2, // plaintext modulus
            p: q * q,
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..1_000 {
            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

            test_mul_relin_opt(&mut rng, &param, m1, m2)?;
        }

        Ok(())
    }

    fn test_mul_relin_opt(mut rng: impl Rng, param: &Param, m1: Rq, m2: Rq) -> Result<()> {
        let (sk, pk) = BFV::new_key(&mut rng, &param)?;

        let rlk = BFV::rlk_key(&mut rng, &param, &sk)?;

        let c1 = BFV::encrypt(&mut rng, &param, &pk, &m1)?;
        let c2 = BFV::encrypt(&mut rng, &param, &pk, &m2)?;

        let c3 = RLWE::mul(param.t, &rlk, &c1, &c2); // uses relinearize internally

        let m3 = BFV::decrypt(&param, &sk, &c3);

        let naive = (m1.clone().to_r() * m2.clone().to_r()).to_rq(param.t); // TODO rm clones
        assert_eq!(
            m3.coeffs().to_vec(),
            naive.coeffs().to_vec(),
            "\n\nfor testing:\nlet m1 = Rq::<T, N>::from_vec_u64(vec!{:?});\nlet m2 = Rq::<T, N>::from_vec_u64(vec!{:?});\n",
            m1.coeffs(),
            m2.coeffs()
            );

        Ok(())
    }
}
