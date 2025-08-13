//! Generalized LWE.
//!

use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, RingParam, Rq, Zq, TR};

use crate::glev::GLev;

// const ERR_SIGMA: f64 = 3.2;
const ERR_SIGMA: f64 = 0.0; // TODO WIP

#[derive(Clone, Copy, Debug)]
pub struct Param {
    pub ring: RingParam,
    pub k: usize,
    pub t: u64,
}
impl Param {
    /// returns the plaintext param
    pub fn pt(&self) -> RingParam {
        // TODO think if maybe return a new truct "PtParam" to differenciate
        // between the ciphertexxt (RingParam) and the plaintext param. Maybe it
        // can be just a wrapper on top of RingParam.
        RingParam {
            q: self.t,
            n: self.ring.n,
        }
    }
    /// returns the LWE param for the given GLWE (self), that is, it uses k=K*N
    /// as the length for the secret key. This follows [2018-421] where
    ///   TLWE sk:  s \in B^n , where n=K*N
    ///   TRLWE sk: s \in B_N[X]^K
    pub fn lwe(&self) -> Self {
        Self {
            ring: RingParam {
                q: self.ring.q,
                n: 1,
            },
            k: self.k * self.ring.n,
            t: self.t,
        }
    }
}

/// GLWE implemented over the `Ring` trait, so that it can be also instantiated
/// over the Torus polynomials 𝕋_<N,q>[X] = 𝕋_q[X]/ (X^N+1).
#[derive(Clone, Debug)]
pub struct GLWE<R: Ring>(pub TR<R>, pub R);

#[derive(Clone, Debug)]
pub struct SecretKey<R: Ring>(pub TR<R>);
#[derive(Clone, Debug)]
pub struct PublicKey<R: Ring>(pub R, pub TR<R>);

// K GLevs, each KSK_i=l GLWEs
#[derive(Clone, Debug)]
pub struct KSK<R: Ring>(Vec<GLev<R>>);

impl<R: Ring> GLWE<R> {
    pub fn zero(k: usize, param: &RingParam) -> Self {
        Self(TR::zero(k, &param), R::zero(&param))
    }
    pub fn from_plaintext(k: usize, param: &RingParam, p: R) -> Self {
        Self(TR::zero(k, &param), p)
    }

    pub fn new_key(mut rng: impl Rng, param: &Param) -> Result<(SecretKey<R>, PublicKey<R>)> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let s: TR<R> = TR::rand(&mut rng, Xi_key, param.k, &param.ring);
        let a: TR<R> = TR::rand(
            &mut rng,
            Uniform::new(0_f64, param.ring.q as f64),
            param.k,
            &param.ring,
        );
        let e = R::rand(&mut rng, Xi_err, &param.ring);

        let pk: PublicKey<R> = PublicKey((&a * &s) + e, a);
        Ok((SecretKey(s), pk))
    }
    pub fn pk_from_sk(mut rng: impl Rng, param: &Param, sk: SecretKey<R>) -> Result<PublicKey<R>> {
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<R> = TR::rand(
            &mut rng,
            Uniform::new(0_f64, param.ring.q as f64),
            param.k,
            &param.ring,
        );
        let e = R::rand(&mut rng, Xi_err, &param.ring);

        let pk: PublicKey<R> = PublicKey((&a * &sk.0) + e, a);
        Ok(pk)
    }

    pub fn new_ksk(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        sk: &SecretKey<R>,
        new_sk: &SecretKey<R>,
    ) -> Result<KSK<R>> {
        debug_assert_eq!(param.k, sk.0.k);
        let k = sk.0.k;
        let r: Vec<GLev<R>> = (0..k)
            .into_iter()
            .map(|i|
                // treat sk_i as the msg being encrypted
                GLev::<R>::encrypt_s(&mut rng, param, beta, l, &new_sk, &sk.0 .r[i]))
            .collect::<Result<Vec<_>>>()?;

        Ok(KSK(r))
    }
    pub fn key_switch(&self, param: &Param, beta: u32, l: u32, ksk: &KSK<R>) -> Self {
        let (a, b): (TR<R>, R) = (self.0.clone(), self.1.clone()); // TODO rm clones

        let lhs: GLWE<R> = GLWE(TR::zero(param.k, &param.ring), b);

        // K iterations, ksk.0 contains K times GLev
        let rhs: GLWE<R> = zip_eq(a.r, ksk.0.clone())
            .map(|(a_i, ksk_i)| ksk_i * a_i.decompose(beta, l)) // dot_product
            .sum();

        lhs - rhs
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(
        mut rng: impl Rng,
        param: &Param,
        sk: &SecretKey<R>,
        m: &R, // already scaled
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<R> = TR::rand(&mut rng, Xi_key, param.k, &param.ring);
        let e = R::rand(&mut rng, Xi_err, &param.ring);

        let b: R = (&a * &sk.0) + m.clone() + e; // TODO rm clone
        Ok(Self(a, b))
    }
    pub fn encrypt(
        mut rng: impl Rng,
        param: &Param,
        pk: &PublicKey<R>,
        m: &R, // already scaled
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u: R = R::rand(&mut rng, Xi_key, &param.ring);

        let e0 = R::rand(&mut rng, Xi_err, &param.ring);
        let e1 = TR::<R>::rand(&mut rng, Xi_err, param.k, &param.ring);

        let b: R = pk.0.clone() * u.clone() + m.clone() + e0; // TODO rm clones
        let d: TR<R> = &pk.1 * &u + e1;

        Ok(Self(d, b))
    }
    // returns m' not downscaled
    pub fn decrypt(&self, sk: &SecretKey<R>) -> R {
        let (d, b): (TR<R>, R) = (self.0.clone(), self.1.clone());
        let p: R = b - &d * &sk.0;
        p
    }
}

// Methods for when Ring=Rq<Q,N>
impl GLWE<Rq> {
    // scale up
    pub fn encode(param: &Param, m: &Rq) -> Rq {
        debug_assert_eq!(param.t, m.param.q);
        let m = m.remodule(param.ring.q);
        let delta = param.ring.q / param.t; // floored
        m * delta
    }
    // scale down
    pub fn decode(param: &Param, m: &Rq) -> Rq {
        let r = m.mul_div_round(param.t, param.ring.q);
        let r: Rq = r.remodule(param.t);
        r
    }
    pub fn mod_switch(&self, p: u64) -> GLWE<Rq> {
        let a: TR<Rq> = TR {
            k: self.0.k,
            r: self.0.r.iter().map(|r| r.mod_switch(p)).collect::<Vec<_>>(),
        };
        let b: Rq = self.1.mod_switch(p);
        GLWE(a, b)
    }
}

impl<R: Ring> Add<GLWE<R>> for GLWE<R> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        debug_assert_eq!(self.0.k, other.0.k);
        debug_assert_eq!(self.1.param(), other.1.param());

        let a: TR<R> = self.0 + other.0;
        let b: R = self.1 + other.1;
        Self(a, b)
    }
}

impl<R: Ring> Add<R> for GLWE<R> {
    type Output = Self;
    fn add(self, plaintext: R) -> Self {
        debug_assert_eq!(self.1.param(), plaintext.param());

        let a: TR<R> = self.0;
        let b: R = self.1 + plaintext;
        Self(a, b)
    }
}
impl<R: Ring> AddAssign for GLWE<R> {
    fn add_assign(&mut self, rhs: Self) {
        debug_assert_eq!(self.0.k, rhs.0.k);
        debug_assert_eq!(self.1.param(), rhs.1.param());

        let k = self.0.k;
        for i in 0..k {
            self.0.r[i] = self.0.r[i].clone() + rhs.0.r[i].clone();
        }
        self.1 = self.1.clone() + rhs.1.clone();
    }
}
impl<R: Ring> Sum<GLWE<R>> for GLWE<R> {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        // let mut acc = GLWE::<R>::zero();
        // for e in iter {
        //     acc += e;
        // }
        // acc
        let first = iter.next().unwrap();
        iter.fold(first, |acc, e| acc + e)
    }
}

impl<R: Ring> Sub<GLWE<R>> for GLWE<R> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        debug_assert_eq!(self.0.k, other.0.k);
        debug_assert_eq!(self.1.param(), other.1.param());

        let a: TR<R> = self.0 - other.0;
        let b: R = self.1 - other.1;
        Self(a, b)
    }
}

impl<R: Ring> Mul<R> for GLWE<R> {
    type Output = Self;
    fn mul(self, plaintext: R) -> Self {
        debug_assert_eq!(self.1.param(), plaintext.param());

        let a: TR<R> = TR {
            k: self.0.k,
            r: self
                .0
                .r
                .iter()
                .map(|r_i| r_i.clone() * plaintext.clone())
                .collect(),
        };
        let b: R = self.1 * plaintext;
        Self(a, b)
    }
}
// for when R = Rq<Q,N>
// impl<const Q: u64, const N: usize, const K: usize> Mul<Rq<Q, N>> for GLWE<Rq<Q, N>, K> {
//     type Output = Self;
//     fn mul(self, plaintext: Rq<Q, N>) -> Self {
//         // first compute the NTT for plaintext, to avoid computing it at each
//         // iteration, speeding up the multiplications
//         let mut plaintext = plaintext.clone();
//         plaintext.compute_evals();
//
//         let a: TR<Rq<Q, N>, K> = TR(self.0 .0.iter().map(|r_i| *r_i * plaintext).collect());
//         let b: Rq<Q, N> = self.1 * plaintext;
//         Self(a, b)
//     }
// }

// impl<R: Ring, const K: usize> Mul<R::C> for GLWE<R, K>
// // where
// //     // R: std::ops::Mul<<R as arith::Ring>::C>,
// //     // Vec<R>: FromIterator<<R as Mul<<R as arith::Ring>::C>>::Output>,
// //     Vec<R>: FromIterator<<R as Mul<<R as arith::Ring>::C>>::Output>,
// {
//     type Output = Self;
//     fn mul(self, e: R::C) -> Self {
//         let a: TR<R, K> = TR(self.0 .0.iter().map(|r_i| *r_i * e.clone()).collect());
//         let b: R = self.1 * e.clone();
//         Self(a, b)
//     }
// }

// impl<const Q: u64, const N: usize, const K: usize> Mul<Zq<Q>> for GLWE<Q, N, K> {
//     type Output = Self;
//     fn mul(self, e: Zq<Q>) -> Self {
//         let a: TR<Rq<Q, N>, K> = TR(self.0 .0.iter().map(|r_i| *r_i * e).collect());
//         let b: Rq<Q, N> = self.1 * e;
//         Self(a, b)
//     }
// }

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt_ring_nq() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 128,
            },
            k: 16,
            t: 32, // plaintext modulus
        };
        // let k: usize = 16;
        type S = GLWE<Rq>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?; // msg
                                                                    // let m: Rq<Q, N> = m.remodule::<Q>();

            let p = S::encode(&param, &m); // plaintext
            let c = S::encrypt(&mut rng, &param, &pk, &p)?; // ciphertext
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode(&param, &p_recovered);

            assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &param, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode(&param, &p_recovered);

            assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));
        }

        Ok(())
    }

    use arith::{Tn, T64};
    pub fn t_encode(param: &RingParam, m: &Rq) -> Tn {
        let p = m.param.q; // plaintext space
        let delta = u64::MAX / p; // floored
        let coeffs = m.coeffs();
        // Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
        // Tn{param, coeffs: array::from_fn(|i| T64(coeffs[i].0 * delta)))
        Tn {
            param: *param,
            coeffs: coeffs.iter().map(|c_i| T64(c_i.v * delta)).collect(),
        }
    }
    pub fn t_decode(param: &Param, pt: &Tn) -> Rq {
        let pt = pt.mul_div_round(param.t, u64::MAX);
        Rq::from_vec_u64(&param.pt(), pt.coeffs().iter().map(|c| c.0).collect())
    }
    #[test]
    fn test_encrypt_decrypt_torus() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: u64::MAX,
                n: 128,
            },
            k: 16,
            t: 32, // plaintext modulus
        };
        type S = GLWE<Tn>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_f64, param.t as f64);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m = Rq::rand(&mut rng, msg_dist, &param.pt()); // msg

            let p = t_encode(&param.ring, &m); // plaintext
            let c = S::encrypt(&mut rng, &param, &pk, &p)?; // ciphertext
            let p_recovered = c.decrypt(&sk);
            let m_recovered = t_decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &param, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = t_decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 128,
            },
            k: 16,
            t: 20, // plaintext modulus
        };
        type S = GLWE<Rq>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Rq = S::encode(&param, &m1); // plaintext
            let p2: Rq = S::encode(&param, &m2); // plaintext

            let c1 = S::encrypt(&mut rng, &param, &pk, &p1)?;
            let c2 = S::encrypt(&mut rng, &param, &pk, &p2)?;

            let c3 = c1 + c2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode(&param, &p3_recovered);

            assert_eq!((m1 + m2).remodule(param.t), m3_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_add_plaintext() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 128,
            },
            k: 16,
            t: 32, // plaintext modulus
        };
        type S = GLWE<Rq>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Rq = S::encode(&param, &m1); // plaintext
            let p2: Rq = S::encode(&param, &m2); // plaintext

            let c1 = S::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode(&param, &p3_recovered);

            assert_eq!((m1 + m2).remodule(param.t), m3_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 16,
            },
            k: 16,
            t: 4, // plaintext modulus
        };
        type S = GLWE<Rq>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Rq = S::encode(&param, &m1); // plaintext
            let p2 = m2.remodule(param.ring.q); // notice we don't encode (scale by delta)

            let c1 = S::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: Rq = c3.decrypt(&sk);
            let m3_recovered: Rq = S::decode(&param, &p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mod_switch() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 8,
            },
            k: 16,
            t: 4, // plaintext modulus, must be a prime or power of a prime
        };
        let new_q: u64 = 2u64.pow(8) + 1;
        // note: wip, Q and P chosen so that P/Q is an integer
        type S = GLWE<Rq>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng, &param)?;

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

            let p = S::encode(&param, &m);
            let c = S::encrypt(&mut rng, &param, &pk, &p)?;

            let c2: GLWE<Rq> = c.mod_switch(new_q);
            assert_eq!(c2.1.param.q, new_q);
            let sk2: SecretKey<Rq> = SecretKey(TR {
                k: param.k,
                r: sk.0.r.iter().map(|s_i| s_i.remodule(new_q)).collect(),
            });

            let p_recovered = c2.decrypt(&sk2);
            let new_param = Param {
                ring: RingParam {
                    q: new_q,
                    n: param.ring.n,
                },
                k: param.k,
                t: param.t,
            };
            let m_recovered = GLWE::<Rq>::decode(&new_param, &p_recovered);

            assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_key_switch() -> Result<()> {
        let param = Param {
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 128,
            },
            k: 16,
            t: 2,
        };
        type S = GLWE<Rq>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng, &param)?;
        let (sk2, _) = S::new_key(&mut rng, &param)?;
        // ksk to switch from sk to sk2
        let ksk = S::new_ksk(&mut rng, &param, beta, l, &sk, &sk2)?;

        let msg_dist = Uniform::new(0_u64, param.t);
        let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
        let p = S::encode(&param, &m); // plaintext
                                       //
        let c = S::encrypt_s(&mut rng, &param, &sk, &p)?;

        let c2 = c.key_switch(&param, beta, l, &ksk);

        // decrypt with the 2nd secret key
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = S::decode(&param, &p_recovered);
        assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));

        // do the same but now encrypting with pk
        let c = S::encrypt(&mut rng, &param, &pk, &p)?;
        let c2 = c.key_switch(&param, beta, l, &ksk);
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = S::decode(&param, &p_recovered);
        assert_eq!(m, m_recovered);

        Ok(())
    }
}
