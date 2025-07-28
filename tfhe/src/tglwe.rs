use anyhow::Result;
use itertools::zip_eq;
use rand::distributions::Standard;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::array;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, Rq, Tn, T64, TR};
use gfhe::{glwe, GLWE};

use crate::tlev::TLev;

pub type SecretKey<const N: usize, const K: usize> = glwe::SecretKey<Tn<N>, K>;
pub type PublicKey<const N: usize, const K: usize> = glwe::PublicKey<Tn<N>, K>;

#[derive(Clone, Debug)]
pub struct TGLWE<const N: usize, const K: usize>(pub GLWE<Tn<N>, K>);

impl<const N: usize, const K: usize> TGLWE<N, K> {
    pub fn zero() -> Self {
        Self(GLWE::<Tn<N>, K>::zero())
    }

    pub fn new_key(rng: impl Rng) -> Result<(SecretKey<N, K>, PublicKey<N, K>)> {
        let (sk, pk) = GLWE::new_key(rng)?;
        // Ok((SecretKey(sk), PublicKey(pk)))
        Ok((sk, pk))
    }

    pub fn encode<const P: u64>(m: &Rq<P, N>) -> Tn<N> {
        let delta = u64::MAX / P; // floored
        let coeffs = m.coeffs();
        Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
    }
    pub fn decode<const P: u64>(p: &Tn<N>) -> Rq<P, N> {
        let p = p.mul_div_round(P, u64::MAX);
        Rq::<P, N>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(rng: impl Rng, sk: &SecretKey<N, K>, p: &Tn<N>) -> Result<Self> {
        let glwe = GLWE::encrypt_s(rng, &sk, p)?;
        Ok(Self(glwe))
    }
    pub fn encrypt(rng: impl Rng, pk: &PublicKey<N, K>, p: &Tn<N>) -> Result<Self> {
        let glwe = GLWE::encrypt(rng, &pk, p)?;
        Ok(Self(glwe))
    }
    pub fn decrypt(&self, sk: &SecretKey<N, K>) -> Tn<N> {
        self.0.decrypt(&sk)
    }
}

impl<const N: usize, const K: usize> Add<TGLWE<N, K>> for TGLWE<N, K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}
impl<const N: usize, const K: usize> AddAssign for TGLWE<N, K> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}
impl<const N: usize, const K: usize> Sum<TGLWE<N, K>> for TGLWE<N, K> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = TGLWE::<N, K>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<const N: usize, const K: usize> Sub<TGLWE<N, K>> for TGLWE<N, K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

// plaintext addition
impl<const N: usize, const K: usize> Add<Tn<N>> for TGLWE<N, K> {
    type Output = Self;
    fn add(self, plaintext: Tn<N>) -> Self {
        let a: TR<Tn<N>, K> = self.0 .0;
        let b: Tn<N> = self.0 .1 + plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext substraction
impl<const N: usize, const K: usize> Sub<Tn<N>> for TGLWE<N, K> {
    type Output = Self;
    fn sub(self, plaintext: Tn<N>) -> Self {
        let a: TR<Tn<N>, K> = self.0 .0;
        let b: Tn<N> = self.0 .1 - plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext multiplication
impl<const N: usize, const K: usize> Mul<Tn<N>> for TGLWE<N, K> {
    type Output = Self;
    fn mul(self, plaintext: Tn<N>) -> Self {
        let a: TR<Tn<N>, K> = TR(self.0 .0 .0.iter().map(|r_i| *r_i * plaintext).collect());
        let b: Tn<N> = self.0 .1 * plaintext;
        Self(GLWE(a, b))
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        const T: u64 = 128; // msg space (msg modulus)
        const N: usize = 64;
        const K: usize = 16;
        type S = TGLWE<N, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p: Tn<N> = S::encode::<T>(&m);

            let c = S::encrypt(&mut rng, &pk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        const T: u64 = 128;
        const N: usize = 64;
        const K: usize = 16;
        type S = TGLWE<N, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<N> = S::encode::<T>(&m1); // plaintext
            let p2: Tn<N> = S::encode::<T>(&m2); // plaintext

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
        const T: u64 = 128;
        const N: usize = 64;
        const K: usize = 16;
        type S = TGLWE<N, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<N> = S::encode::<T>(&m1); // plaintext
            let p2: Tn<N> = S::encode::<T>(&m2); // plaintext

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);

            assert_eq!(m1 + m2, m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        const T: u64 = 128;
        const N: usize = 64;
        const K: usize = 16;
        type S = TGLWE<N, K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<N> = S::encode::<T>(&m1);
            // don't scale up p2, set it directly from m2
            let p2: Tn<N> = Tn(array::from_fn(|i| T64(m2.coeffs()[i].0)));

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: Tn<N> = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), m3_recovered);
        }

        Ok(())
    }
}
