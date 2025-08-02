//! Generalized LWE.
//!

use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, Rq, Zq, TR};

use crate::glev::GLev;

// const ERR_SIGMA: f64 = 3.2;
const ERR_SIGMA: f64 = 0.0; // TODO WIP

/// GLWE implemented over the `Ring` trait, so that it can be also instantiated
/// over the Torus polynomials 𝕋_<N,q>[X] = 𝕋_q[X]/ (X^N+1).
#[derive(Clone, Debug)]
pub struct GLWE<R: Ring, const K: usize>(pub TR<R, K>, pub R);

#[derive(Clone, Debug)]
pub struct SecretKey<R: Ring, const K: usize>(pub TR<R, K>);
#[derive(Clone, Debug)]
pub struct PublicKey<R: Ring, const K: usize>(pub R, pub TR<R, K>);

// K GLevs, each KSK_i=l GLWEs
#[derive(Clone, Debug)]
pub struct KSK<R: Ring, const K: usize>(Vec<GLev<R, K>>);

impl<R: Ring, const K: usize> GLWE<R, K> {
    pub fn zero() -> Self {
        Self(TR::zero(), R::zero())
    }

    pub fn new_key(mut rng: impl Rng) -> Result<(SecretKey<R, K>, PublicKey<R, K>)> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let s: TR<R, K> = TR::rand(&mut rng, Xi_key);
        let a: TR<R, K> = TR::rand(&mut rng, Uniform::new(0_f64, R::Q as f64));
        let e = R::rand(&mut rng, Xi_err);

        let pk: PublicKey<R, K> = PublicKey((&a * &s) + e, a);
        Ok((SecretKey(s), pk))
    }
    pub fn pk_from_sk(mut rng: impl Rng, sk: SecretKey<R, K>) -> Result<PublicKey<R, K>> {
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<R, K> = TR::rand(&mut rng, Uniform::new(0_f64, R::Q as f64));
        let e = R::rand(&mut rng, Xi_err);

        let pk: PublicKey<R, K> = PublicKey((&a * &sk.0) + e, a);
        Ok(pk)
    }

    pub fn new_ksk(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<R, K>,
        new_sk: &SecretKey<R, K>,
    ) -> Result<KSK<R, K>> {
        let r: Vec<GLev<R, K>> = (0..K)
            .into_iter()
            .map(|i|
                // treat sk_i as the msg being encrypted
                GLev::<R, K>::encrypt_s(&mut rng, beta, l, &new_sk, &sk.0 .0[i]))
            .collect::<Result<Vec<_>>>()?;

        Ok(KSK(r))
    }
    pub fn key_switch(&self, beta: u32, l: u32, ksk: &KSK<R, K>) -> Self {
        let (a, b): (TR<R, K>, R) = (self.0.clone(), self.1);

        let lhs: GLWE<R, K> = GLWE(TR::zero(), b);

        // K iterations, ksk.0 contains K times GLev
        let rhs: GLWE<R, K> = zip_eq(a.0, ksk.0.clone())
            .map(|(a_i, ksk_i)| ksk_i * a_i.decompose(beta, l)) // dot_product
            .sum();

        lhs - rhs
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(
        mut rng: impl Rng,
        sk: &SecretKey<R, K>,
        m: &R, // already scaled
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<R, K> = TR::rand(&mut rng, Xi_key);
        let e = R::rand(&mut rng, Xi_err);

        let b: R = (&a * &sk.0) + *m + e;
        Ok(Self(a, b))
    }
    pub fn encrypt(
        mut rng: impl Rng,
        pk: &PublicKey<R, K>,
        m: &R, // already scaled
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u: R = R::rand(&mut rng, Xi_key);

        let e0 = R::rand(&mut rng, Xi_err);
        let e1 = TR::<R, K>::rand(&mut rng, Xi_err);

        let b: R = pk.0.clone() * u.clone() + *m + e0;
        let d: TR<R, K> = &pk.1 * &u + e1;

        Ok(Self(d, b))
    }
    // returns m' not downscaled
    pub fn decrypt(&self, sk: &SecretKey<R, K>) -> R {
        let (d, b): (TR<R, K>, R) = (self.0.clone(), self.1);
        let p: R = b - &d * &sk.0;
        p
    }
}

// Methods for when Ring=Rq<Q,N>
impl<const Q: u64, const N: usize, const K: usize> GLWE<Rq<Q, N>, K> {
    // scale up
    pub fn encode<const T: u64>(m: &Rq<T, N>) -> Rq<Q, N> {
        let m = m.remodule::<Q>();
        let delta = Q / T; // floored
        m * delta
    }
    // scale down
    pub fn decode<const T: u64>(m: &Rq<Q, N>) -> Rq<T, N> {
        let r = m.mul_div_round(T, Q);
        let r: Rq<T, N> = r.remodule::<T>();
        r
    }
    pub fn mod_switch<const P: u64>(&self) -> GLWE<Rq<P, N>, K> {
        let a: TR<Rq<P, N>, K> = TR(self
            .0
             .0
            .iter()
            .map(|r| r.mod_switch::<P>())
            .collect::<Vec<_>>());
        let b: Rq<P, N> = self.1.mod_switch::<P>();
        GLWE(a, b)
    }
}

impl<R: Ring, const K: usize> Add<GLWE<R, K>> for GLWE<R, K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let a: TR<R, K> = self.0 + other.0;
        let b: R = self.1 + other.1;
        Self(a, b)
    }
}

impl<R: Ring, const K: usize> Add<R> for GLWE<R, K> {
    type Output = Self;
    fn add(self, plaintext: R) -> Self {
        let a: TR<R, K> = self.0;
        let b: R = self.1 + plaintext;
        Self(a, b)
    }
}
impl<R: Ring, const K: usize> AddAssign for GLWE<R, K> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..K {
            self.0 .0[i] = self.0 .0[i].clone() + rhs.0 .0[i].clone();
        }
        self.1 = self.1.clone() + rhs.1.clone();
    }
}
impl<R: Ring, const K: usize> Sum<GLWE<R, K>> for GLWE<R, K> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = GLWE::<R, K>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<R: Ring, const K: usize> Sub<GLWE<R, K>> for GLWE<R, K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let a: TR<R, K> = self.0 - other.0;
        let b: R = self.1 - other.1;
        Self(a, b)
    }
}

impl<R: Ring, const K: usize> Mul<R> for GLWE<R, K> {
    type Output = Self;
    fn mul(self, plaintext: R) -> Self {
        let a: TR<R, K> = TR(self.0 .0.iter().map(|r_i| *r_i * plaintext).collect());
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
    fn test_encrypt_decrypt() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 32; // plaintext modulus
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?; // msg
                                                               // let m: Rq<Q, N> = m.remodule::<Q>();

            let p = S::encode::<T>(&m); // plaintext
            let c = S::encrypt(&mut rng, &pk, &p)?; // ciphertext
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }

    use arith::{Tn, T64};
    use std::array;
    pub fn t_encode<const P: u64>(m: &Rq<P, 4>) -> Tn<4> {
        let delta = u64::MAX / P; // floored
        let coeffs = m.coeffs();
        Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
    }
    pub fn t_decode<const P: u64>(p: &Tn<4>) -> Rq<P, 4> {
        let p = p.mul_div_round(P, u64::MAX);
        Rq::<P, 4>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }
    #[test]
    fn test_encrypt_decrypt_torus() -> Result<()> {
        const N: usize = 128;
        const T: u64 = 32; // plaintext modulus
        const K: usize = 16;
        type S = GLWE<Tn<4>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_f64, T as f64);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m = Rq::<T, 4>::rand(&mut rng, msg_dist); // msg

            let p = t_encode::<T>(&m); // plaintext
            let c = S::encrypt(&mut rng, &pk, &p)?; // ciphertext
            let p_recovered = c.decrypt(&sk);
            let m_recovered = t_decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = t_decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 20;
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Rq<Q, N> = S::encode::<T>(&m1); // plaintext
            let p2: Rq<Q, N> = S::encode::<T>(&m2); // plaintext

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;
            let c2 = S::encrypt(&mut rng, &pk, &p2)?;

            let c3 = c1 + c2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);

            assert_eq!((m1 + m2).remodule::<T>(), m3_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_add_plaintext() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 32;
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Rq<Q, N> = S::encode::<T>(&m1); // plaintext
            let p2: Rq<Q, N> = S::encode::<T>(&m2); // plaintext

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);

            assert_eq!((m1 + m2).remodule::<T>(), m3_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 16;
        const T: u64 = 4;
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Rq<Q, N> = S::encode::<T>(&m1); // plaintext
            let p2 = m2.remodule::<Q>(); // notice we don't encode (scale by delta)

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: Rq<Q, N> = c3.decrypt(&sk);
            let m3_recovered: Rq<T, N> = S::decode::<T>(&p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mod_switch() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const P: u64 = 2u64.pow(8) + 1;
        // note: wip, Q and P chosen so that P/Q is an integer
        const N: usize = 8;
        const T: u64 = 4; // plaintext modulus, must be a prime or power of a prime
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;

            let p = S::encode::<T>(&m);
            let c = S::encrypt(&mut rng, &pk, &p)?;

            let c2: GLWE<Rq<P, N>, K> = c.mod_switch::<P>();
            let sk2: SecretKey<Rq<P, N>, K> =
                SecretKey(TR(sk.0 .0.iter().map(|s_i| s_i.remodule::<P>()).collect()));

            let p_recovered = c2.decrypt(&sk2);
            let m_recovered = GLWE::<Rq<P, N>, K>::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_key_switch() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 2; // plaintext modulus
        const K: usize = 16;
        type S = GLWE<Rq<Q, N>, K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;
        let (sk2, _) = S::new_key(&mut rng)?;
        // ksk to switch from sk to sk2
        let ksk = S::new_ksk(&mut rng, beta, l, &sk, &sk2)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let p = S::encode::<T>(&m); // plaintext
                                    //
        let c = S::encrypt_s(&mut rng, &sk, &p)?;

        let c2 = c.key_switch(beta, l, &ksk);

        // decrypt with the 2nd secret key
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = S::decode::<T>(&p_recovered);
        assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());

        // do the same but now encrypting with pk
        let c = S::encrypt(&mut rng, &pk, &p)?;
        let c2 = c.key_switch(beta, l, &ksk);
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = S::decode::<T>(&p_recovered);
        assert_eq!(m, m_recovered);

        Ok(())
    }
}
