//! Implementation of BFV https://eprint.iacr.org/2012/144.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

use anyhow::{anyhow, Result};
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::ops;

use arithmetic::{Rq, Zq, R};

// error deviation for the Gaussian(Normal) distribution
// sigma=3.2 from: https://eprint.iacr.org/2022/162.pdf page 5
const ERR_SIGMA: f64 = 3.2;
// const ERR_SIGMA: f64 = 0.0; // TODO WIP

#[derive(Clone, Debug)]
pub struct SecretKey<const Q: u64, const N: usize>(Rq<Q, N>);

#[derive(Clone, Debug)]
pub struct PublicKey<const Q: u64, const N: usize>(Rq<Q, N>, Rq<Q, N>);

/// Relinearization key
#[derive(Clone, Debug)]
pub struct RLK<const PQ: u64, const N: usize>(Rq<PQ, N>, Rq<PQ, N>);

// impl<const PQ: u64, const N: usize> RLK<Q, PQ, N> {
//     // const P: u64 = PQ / Q;
//
//     // const PQ: u64 = P * Q;
// }

// RLWE ciphertext
#[derive(Clone, Debug)]
pub struct RLWE<const Q: u64, const N: usize>(Rq<Q, N>, Rq<Q, N>);

impl<const Q: u64, const N: usize> RLWE<Q, N> {
    fn add(lhs: Self, rhs: Self) -> Self {
        RLWE::<Q, N>(lhs.0 + rhs.0, lhs.1 + rhs.1)
    }
    pub fn remodule<const P: u64>(&self) -> RLWE<P, N> {
        let x = self.0.remodule::<P>();
        let y = self.1.remodule::<P>();
        RLWE::<P, N>(x, y)
    }
    fn tensor<const PQ: u64, const T: u64>(a: &Self, b: &Self) -> (Rq<Q, N>, Rq<Q, N>, Rq<Q, N>) {
        // expand Q->PQ // TODO rm
        // get the coefficients in Z, ie. interpret a,b \in R (instead of R_q)
        let a0: R<N> = a.0.to_r();
        let a1: R<N> = a.1.to_r();
        let b0: R<N> = b.0.to_r();
        let b1: R<N> = b.1.to_r();

        // tensor (\in R)
        use arithmetic::ring::naive_mul;
        let c0: Vec<i64> = naive_mul(&a0, &b0);
        let c1_l: Vec<i64> = naive_mul(&a0, &b1);
        let c1_r = naive_mul(&a1, &b0);
        let c1: Vec<i64> = itertools::zip_eq(c1_l, c1_r).map(|(l, r)| l + r).collect();
        let c2: Vec<i64> = naive_mul(&a1, &b1);

        // scale down, module Q, so result is \in R_q
        let c0: Rq<Q, N> = arithmetic::ring::mul_div_round::<Q, N>(c0, T, Q);
        let c1: Rq<Q, N> = arithmetic::ring::mul_div_round::<Q, N>(c1, T, Q);
        let c2: Rq<Q, N> = arithmetic::ring::mul_div_round::<Q, N>(c2, T, Q);

        (c0, c1, c2)
    }
    fn tensor_DBG<const PQ: u64, const T: u64>(
        a: &Self,
        b: &Self,
    ) -> (Rq<Q, N>, Rq<Q, N>, Rq<Q, N>) {
        // iacr 2021-204:
        // expand Q->PQ
        // let a: RLWE<PQ, N> = a.remodule::<PQ>();
        // let b: RLWE<PQ, N> = b.remodule::<PQ>();

        // tensor
        let c0: Rq<Q, N> = a.0 * b.0; // NTT mul
        let c1: Rq<Q, N> = a.0 * b.1 + a.1 * b.0; // NTT mul
        let c2: Rq<Q, N> = a.1 * b.1; // NTT mul

        // scale down
        let c0: Rq<Q, N> = c0.mul_div_round(T, Q);
        let c1: Rq<Q, N> = c1.mul_div_round(T, Q);
        let c2: Rq<Q, N> = c2.mul_div_round(T, Q);

        // expand^-1 PQ->Q
        // let c0: Rq<Q, N> = c0.remodule::<Q>();
        // let c1: Rq<Q, N> = c1.remodule::<Q>();
        // let c2: Rq<Q, N> = c2.remodule::<Q>();

        (c0, c1, c2)
    }
    /// ciphertext multiplication
    fn mul<const PQ: u64, const T: u64>(rlk: &RLK<PQ, N>, a: &Self, b: &Self) -> Self {
        let (c0, c1, c2) = Self::tensor::<PQ, T>(a, b);
        BFV::<Q, N, T>::relinearize_204::<PQ>(&rlk, &c0, &c1, &c2)
    }
}
// naive mul in the ring Rq, reusing the ring::naive_mul and then applying mod(X^N +1)
fn tmp_naive_mul<const Q: u64, const N: usize>(a: Rq<Q, N>, b: Rq<Q, N>) -> Rq<Q, N> {
    Rq::<Q, N>::from_vec_i64(arithmetic::ring::naive_mul(&a.to_r(), &b.to_r()))
}

impl<const Q: u64, const N: usize> ops::Add<RLWE<Q, N>> for RLWE<Q, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::add(self, rhs)
    }
}

impl<const Q: u64, const N: usize, const T: u64> ops::Add<&Rq<T, N>> for &RLWE<Q, N> {
    type Output = RLWE<Q, N>;
    fn add(self, rhs: &Rq<T, N>) -> Self::Output {
        BFV::<Q, N, T>::add_const(self, rhs)
    }
}
impl<const Q: u64, const N: usize, const T: u64> ops::Mul<&Rq<T, N>> for &RLWE<Q, N> {
    type Output = RLWE<Q, N>;
    fn mul(self, rhs: &Rq<T, N>) -> Self::Output {
        BFV::<Q, N, T>::mul_const(&self, rhs)
    }
}

pub struct BFV<const Q: u64, const N: usize, const T: u64> {}

impl<const Q: u64, const N: usize, const T: u64> BFV<Q, N, T> {
    const DELTA: u64 = Q / T;

    /// generate a new key pair (privK, pubK)
    pub fn new_key(mut rng: impl Rng) -> Result<(SecretKey<Q, N>, PublicKey<Q, N>)> {
        // WIP: review probabilities

        // let Xi_key = Uniform::new(-1_f64, 1_f64);
        let Xi_key = Uniform::new(0_u64, 2_u64);
        // let Xi_key = Uniform::new(0_u64, 2_u64);
        // use rand::distributions::Bernoulli;
        // let Xi_key = Bernoulli::new(0.5)?;
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        // secret key
        let s = Rq::<Q, N>::rand_u64(&mut rng, Xi_key)?;

        #[cfg(test)] // sanity check
        assert!(s.infinity_norm() <= 1, "{:?}", s.coeffs());

        // pk = (-a * s + e, a)
        let a = Rq::<Q, N>::rand_u64(&mut rng, Uniform::new(0_u64, Q))?;
        let e = Rq::<Q, N>::rand_f64(&mut rng, Xi_err)?;
        let pk: PublicKey<Q, N> = PublicKey((&(-a) * &s) + e, a.clone());
        Ok((SecretKey(s), pk))
    }

    pub fn encrypt(mut rng: impl Rng, pk: &PublicKey<Q, N>, m: &Rq<T, N>) -> Result<RLWE<Q, N>> {
        // let Xi_key = Uniform::new(-1_f64, 1_f64);
        // let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_key = Uniform::new(0_u64, 2_u64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u = Rq::<Q, N>::rand_u64(&mut rng, Xi_key)?;
        let e_1 = Rq::<Q, N>::rand_f64_abs(&mut rng, Xi_err)?;
        let e_2 = Rq::<Q, N>::rand_f64_abs(&mut rng, Xi_err)?;

        // println!("{:?}", &e_1.coeffs());
        // println!("{:?}", &e_2.coeffs());

        #[cfg(test)] // sanity check
        assert!(u.infinity_norm() <= 1, "{:?}", u.coeffs());

        // migrate m's coeffs to the bigger modulus Q (from T)
        let m = m.remodule::<Q>();
        #[cfg(test)]
        {
            // sanity check // TODO rm
            let m_remod_naive =
                Rq::<Q, N>::from_vec_u64(m.coeffs().iter().map(|m_i| m_i.0).collect());
            assert_eq!(m_remod_naive, m);
        }

        // let c0 = &pk.0 * &u + e_1 + m * Self::DELTA;
        let c0 = &pk.0 * &u + e_1 + m.mul_div_round(Q, T); // TODO use DELTA?
        let c1 = &pk.1 * &u + e_2;
        // let c0 = tmp_naive_mul(pk.0, u) + e_1 + m * Self::DELTA;
        // let c0 = tmp_naive_mul(pk.0, u) + e_1 + m.mul_div_round(Q, T);
        // let c1 = tmp_naive_mul(pk.1, u)
        //     // &pk.1 * &u
        //     + e_2;
        Ok(RLWE::<Q, N>(c0, c1))
    }

    pub fn decrypt(sk: &SecretKey<Q, N>, c: &RLWE<Q, N>) -> Rq<T, N> {
        let cs = c.0 + c.1 * sk.0; // done in mod q

        // let c1s = tmp_naive_mul(c.1, sk.0);
        // // let c1s = arithmetic::ring::naive_mul(&c.1.to_r(), &sk.0.to_r()); // TODO rm
        // // let c1s = Rq::<Q, N>::from_vec_i64(c1s);
        // let cs = c.0 + c1s;

        // let r: Vec<u64> = cs
        //     .coeffs()
        //     .iter()
        //     .map(|e| ((T as f64 * e.0 as f64) / Q as f64).round() as u64)
        //     .collect();
        // Rq::<T, N>::from_vec_u64(r)

        let r: Rq<Q, N> = cs.mul_div_round(T, Q);
        r.remodule::<T>()
    }

    fn add_const(c: &RLWE<Q, N>, m: &Rq<T, N>) -> RLWE<Q, N> {
        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule::<Q>();
        RLWE::<Q, N>(c.0 + m * Self::DELTA, c.1)
    }
    fn mul_const(c: &RLWE<Q, N>, m: &Rq<T, N>) -> RLWE<Q, N> {
        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule::<Q>();
        RLWE::<Q, N>(c.0 * m * Self::DELTA, c.1)
    }

    fn rlk_key<const PQ: u64>(mut rng: impl Rng, s: &SecretKey<Q, N>) -> Result<RLK<PQ, N>> {
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?; // TODO review Xi' instead of Xi

        let s = s.0.remodule::<PQ>();
        let a = Rq::<PQ, N>::rand_u64(&mut rng, Uniform::new(0_u64, PQ))?;
        let e = Rq::<PQ, N>::rand_f64(&mut rng, Xi_err)?;

        // let rlk_1: Rq<PQ, N> = (&(-a) * &s) - e + (s * s) * P;
        let P = PQ / Q;

        // let rlk: RLK<PQ, N> = RLK::<PQ, N>((&(-a) * &s) - e + (s * s) * P, a.clone());
        let rlk: RLK<PQ, N> = RLK::<PQ, N>(
            -(tmp_naive_mul(a, s) + e) + tmp_naive_mul(s, s) * P,
            a.clone(),
        );
        // let rlk: RLK<PQ, N> = RLK::<PQ, N>(-(&a * &s + e) + (s * s) * P, a.clone());

        // let pq = P * Q;
        // let a = Rq::<Q, N>::rand_u64(&mut rng, Uniform::new(0_u64, pq))?;
        // let e = Rq::<Q, N>::rand_f64(&mut rng, Xi_err)?;
        //
        // let rlk_0: Rq<Q, N> = (&(-a) * &s) - e + (s * s) * P;
        // let rlk_0 = rlk_0.remodule::<>();
        // let rlk: RLK<Q, N> = RLK(rlk_0, a);
        Ok(rlk)
    }

    fn relinearize<const PQ: u64>(
        rlk: &RLK<PQ, N>,
        c0: &Rq<Q, N>,
        c1: &Rq<Q, N>,
        c2: &Rq<Q, N>,
    ) -> RLWE<Q, N> {
        let P = PQ / Q;

        // let c2 = c2.remodule::<PQ>();
        // let c2 = c2.to_r();

        let c2rlk0: Vec<f64> = (c2.to_r() * rlk.0.to_r())
            .coeffs()
            .iter()
            .map(|e| (*e as f64 / P as f64).round())
            .collect();

        let c2rlk1: Vec<f64> = (c2.to_r() * rlk.1.to_r())
            .coeffs()
            .iter()
            .map(|e| (*e as f64 / P as f64).round())
            .collect();

        let r0 = Rq::<Q, N>::from_vec_f64(c2rlk0);
        let r1 = Rq::<Q, N>::from_vec_f64(c2rlk1);

        let res = RLWE::<Q, N>(c0 + &r0, c1 + &r1);
        res
    }
    fn relinearize_204<const PQ: u64>(
        rlk: &RLK<PQ, N>,
        c0: &Rq<Q, N>,
        c1: &Rq<Q, N>,
        c2: &Rq<Q, N>,
    ) -> RLWE<Q, N> {
        let P = PQ / Q;

        // let c2 = c2.remodule::<PQ>();
        // let c2 = c2.to_r();

        // let c2rlk0: Vec<f64> = (c2.remodule::<PQ>() * rlk.0)
        use arithmetic::ring::naive_mul;
        let c2rlk0: Vec<i64> = naive_mul(&c2.to_r(), &rlk.0.to_r());
        // .coeffs()
        // .iter()
        // .map(|e| (*e as f64 / P as f64).round())
        // .collect();

        let c2rlk1: Vec<i64> = naive_mul(&c2.to_r(), &rlk.1.to_r());
        // .coeffs()
        // .iter()
        // .map(|e| (*e as f64 / P as f64).round())
        // .collect();
        //
        let r0: Rq<Q, N> = arithmetic::ring::mul_div_round::<Q, N>(c2rlk0, 1, P);
        let r1: Rq<Q, N> = arithmetic::ring::mul_div_round::<Q, N>(c2rlk1, 1, P);

        // let r0 = Rq::<Q, N>::from_vec_f64(c2rlk0);
        // let r1 = Rq::<Q, N>::from_vec_f64(c2rlk1);

        let res = RLWE::<Q, N>(c0 + &r0, c1 + &r1);
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
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 32;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

            let c = S::encrypt(&mut rng, &pk, &m)?;
            let m_recovered = S::decrypt(&sk, &c);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        // const N: usize = 32;
        const N: usize = 4;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

        let c1 = S::encrypt(&mut rng, &pk, &m1)?;
        let c2 = S::encrypt(&mut rng, &pk, &m2)?;

        let c3 = c1 + c2;

        let m3_recovered = S::decrypt(&sk, &c3);

        assert_eq!(m1 + m2, m3_recovered);

        Ok(())
    }

    #[test]
    fn test_constant_add() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 32;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let m2_const = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

        let c1 = S::encrypt(&mut rng, &pk, &m1)?;

        let c3_add = &c1 + &m2_const;
        // let c3_mul = &c1 * &m2_const;

        let m3_add_recovered = S::decrypt(&sk, &c3_add);
        // let m3_mul_recovered = S::decrypt(&sk, &c3_mul);

        assert_eq!(m1 + m2_const, m3_add_recovered);
        //
        // let mut mul_res = naive_poly_mul::<T>(&m1.coeffs().to_vec(), &m2_const.coeffs().to_vec());
        // arithmetic::ring::modulus::<T, N>(&mut mul_res);
        // dbg!(&mul_res);
        // let mul_res_2 =
        //     naive_poly_mul_2::<T, N>(&m1.coeffs().to_vec(), &m2_const.coeffs().to_vec());
        // assert_eq!(mul_res, mul_res_2);
        // let mul_res = PR::<T, N>::from_vec(mul_res);
        // assert_eq!(mul_res.coeffs(), m3_mul_recovered.coeffs());

        Ok(())
    }

    #[test]
    fn test_tensor() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1; // q prime, and 2^q + 1 shape
        const N: usize = 8;
        const T: u64 = 4; // plaintext modulus

        // const P: u64 = Q;
        const P: u64 = Q * Q;
        // const P: u64 = 2_u64.pow(13) * Q + 1;
        const PQ: u64 = P * Q;

        let mut rng = rand::thread_rng();

        let msg_dist = Uniform::new(0_u64, T);
        for _ in 0..10_000 {
            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

            test_tensor_opt::<Q, N, T, PQ>(&mut rng, m1, m2)?;
        }

        Ok(())
    }
    fn test_tensor_opt<const Q: u64, const N: usize, const T: u64, const PQ: u64>(
        mut rng: impl Rng,
        m1: Rq<T, N>,
        m2: Rq<T, N>,
    ) -> Result<()> {
        let (sk, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;

        let c1 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m1)?;
        let c2 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m2)?;

        let (c_a, c_b, c_c) = RLWE::<Q, N>::tensor::<PQ, T>(&c1, &c2);

        // decrypt non-relinearized mul result
        let m3: Rq<Q, N> = c_a + c_b * sk.0 + c_c * sk.0 * sk.0;
        let m3: Rq<Q, N> = m3.mul_div_round(T, Q);
        let m3 = m3.remodule::<T>();

        let naive = (m1.to_r() * m2.to_r()).to_rq::<T>();
        assert_eq!(
            m3.coeffs().to_vec(),
            naive.coeffs().to_vec(),
            "\n\nfor testing:\nlet m1 = Rq::<T, N>::from_vec_u64(vec!{:?});\nlet m2 = Rq::<T, N>::from_vec_u64(vec!{:?});\n",
            m1.coeffs(),
            m2.coeffs()
        );
        if m3.coeffs().to_vec() != naive.coeffs().to_vec() {
            return Err(anyhow!("not eq"));
        }

        Ok(())
    }

    fn test_tensor_opt_DBG<const Q: u64, const N: usize, const T: u64, const PQ: u64>(
        mut rng: impl Rng,
        m1: Rq<T, N>,
        m2: Rq<T, N>,
    ) -> Result<()> {
        let (sk, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;

        let c1 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m1)?;
        let c2 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m2)?;

        let (c_a, c_b, c_c) = RLWE::<Q, N>::tensor::<PQ, T>(&c1, &c2);

        // decrypt non-relinearized mul result
        let m3: Rq<Q, N> = c_a + c_b * sk.0 + c_c * sk.0 * sk.0;
        dbg!(m3);
        let m3: Rq<Q, N> = m3.mul_div_round(T, Q);
        dbg!(m3);
        let m3 = m3.remodule::<T>();
        dbg!(m3);

        // let naive = (m1.to_r() * m2.to_r()).to_rq::<T>();
        // let naive = m1.remodule::<Q>() * m2.remodule::<Q>();
        let naive = (m1.remodule::<Q>() * m2.remodule::<Q>()).remodule::<T>();
        dbg!(naive);
        assert_eq!(
            m3.coeffs().to_vec(),
            naive.coeffs().to_vec(),
            "\n\nfor testing:\nlet m1 = Rq::<T, N>::from_vec_u64(vec!{:?});\nlet m2 = Rq::<T, N>::from_vec_u64(vec!{:?});\n",
            m1.coeffs(),
            m2.coeffs()
        );
        // if m3.coeffs().to_vec() != naive.coeffs().to_vec() {
        //     return Err(anyhow!("not eq"));
        // }

        Ok(())
    }

    #[test]
    fn test_mul_relin() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 32;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        const P: u64 = Q * Q;
        const PQ: u64 = P * Q;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..100 {
            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

            test_mul_relin_opt::<Q, N, T, PQ>(&mut rng, m1, m2)?;
        }

        Ok(())
    }

    fn test_mul_relin_opt<const Q: u64, const N: usize, const T: u64, const PQ: u64>(
        mut rng: impl Rng,
        m1: Rq<T, N>,
        m2: Rq<T, N>,
    ) -> Result<()> {
        let (sk, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;

        let rlk = BFV::<Q, N, T>::rlk_key::<PQ>(&mut rng, &sk)?;

        let c1 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m1)?;
        let c2 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &m2)?;

        let c3 = RLWE::<Q, N>::mul::<PQ, T>(&rlk, &c1, &c2);

        let m3 = BFV::<Q, N, T>::decrypt(&sk, &c3);

        let naive = (m1.to_r() * m2.to_r()).to_rq::<T>();
        assert_eq!(m3.coeffs().to_vec(), naive.coeffs().to_vec(),
            "\n\nfor testing:\nlet m1 = Rq::<T, N>::from_vec_u64(vec!{:?});\nlet m2 = Rq::<T, N>::from_vec_u64(vec!{:?});\n",
            m1.coeffs(),
            m2.coeffs()
            );

        Ok(())
    }

    #[test]
    fn test_naive_mul() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1; // prime, and 2^q + 1 shape
        const N: usize = 4;
        const T: u64 = 4; // plaintext modulus

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);
        // for _ in 0..10_000 {
        for _ in 0..2 {
            println!("---");
            // let a = Rq::<Q, N>::rand_u64(&mut rng, msg_dist)?;
            // let b = Rq::<Q, N>::rand_u64(&mut rng, msg_dist)?;
            // let a = Rq::<Q, N>::from_vec_u64(vec![Q - 1, Q - 2, Q - 3, Q - 3]);
            // let b = Rq::<Q, N>::from_vec_u64(vec![Q - 3, Q - 3, Q - 2, Q - 1]);
            let a = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let b = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            dbg!(&a);
            dbg!(&b);
            let (_, pk) = BFV::<Q, N, T>::new_key(&mut rng)?;
            let ciph1 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &a)?;
            let ciph2 = BFV::<Q, N, T>::encrypt(&mut rng, &pk, &b)?;
            let a = ciph1.0;
            let b = ciph2.0;
            dbg!(&a);
            dbg!(&b);

            let c0: Vec<i64> = arithmetic::ring::naive_mul(&a.to_r(), &b.to_r());
            let c0 = Rq::<Q, N>::from_vec_i64(c0);
            let c1 = tmp_naive_mul(a, b); // naive mul
            let c2: Rq<Q, N> = a * b; // NTT mul

            println!("{:?}", c0.coeffs());
            println!("{:?}", c1.coeffs());
            println!("{:?}", c2.coeffs());
            assert_eq!(c0, c2);
            assert_eq!(c1, c2);

            // scale by Delta=Q/T
            let a = a.mul_div_round(Q, T);
            let b = b.mul_div_round(Q, T);
            dbg!(&a);
            dbg!(&b);

            let c0: Vec<i64> = arithmetic::ring::naive_mul(&a.to_r(), &b.to_r());
            let c0 = Rq::<Q, N>::from_vec_i64(c0);
            let c1 = tmp_naive_mul(a, b); // naive mul
            let c2: Rq<Q, N> = a * b; // NTT mul

            println!("{:?}", c0.coeffs());
            println!("{:?}", c1.coeffs());
            println!("{:?}", c2.coeffs());
            assert_eq!(c0, c2);
            assert_eq!(c1, c2);

            let c0 = c0.mul_div_round(T, Q);
            let c1 = c1.mul_div_round(T, Q);
            let c2 = c2.mul_div_round(T, Q);
            println!("{:?}", c0.coeffs());
            println!("{:?}", c1.coeffs());
            println!("{:?}", c2.coeffs());
            assert_eq!(c0, c2);
            assert_eq!(c1, c2);

            /*
            // now same as before, but multiplying by T/Q
            let c0: Vec<i64> = arithmetic::ring::naive_mul(&a.to_r(), &b.to_r());
            let c0: Vec<f64> = c0
                .iter()
                .map(|e| ((T as f64 * *e as f64) / Q as f64).round())
                .collect();
            let c0 = Rq::<Q, N>::from_vec_f64(c0);
            dbg!(&c0.coeffs());

            let a = a.mul_div_round(T, Q);
            let b = b.mul_div_round(T, Q);
            println!("a{:?}", a.coeffs());
            println!("b{:?}", b.coeffs());

            let c4: Vec<i64> = arithmetic::ring::naive_mul(&a.to_r(), &b.to_r());
            let c4 = Rq::<Q, N>::from_vec_i64(c4);
            let c4 = c4.mul_div_round(T, Q);

            let c1 = tmp_naive_mul(a, b); // naive mul
            let c1 = c1.mul_div_round(T, Q);
            let c2 = a * b; // NTT mul
            let c2 = c2.mul_div_round(T, Q);

            println!("{:?}", c0.coeffs());
            println!("{:?}", c1.coeffs());
            println!("{:?}", c2.coeffs());
            println!("{:?}", c4.coeffs());
            assert_eq!(c0, c2);
            assert_eq!(c4, c2);
            assert_eq!(c1, c2);
             */
        }

        Ok(())
    }
}
