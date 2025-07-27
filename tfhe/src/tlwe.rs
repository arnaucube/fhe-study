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

// #[derive(Clone, Debug)]
// pub struct SecretKey<const K: usize>(glwe::SecretKey<T64, K>);
pub type SecretKey<const K: usize> = glwe::SecretKey<T64, K>;

// #[derive(Clone, Debug)]
// pub struct PublicKey<const K: usize>(glwe::PublicKey<T64, K>);
pub type PublicKey<const K: usize> = glwe::PublicKey<T64, K>;

#[derive(Clone, Debug)]
pub struct TLWE<const K: usize>(pub GLWE<T64, K>);

impl<const K: usize> TLWE<K> {
    pub fn zero() -> Self {
        Self(GLWE::<T64, K>::zero())
    }

    pub fn new_key(rng: impl Rng) -> Result<(SecretKey<K>, PublicKey<K>)> {
        let (sk, pk) = GLWE::new_key(rng)?;
        // Ok((SecretKey(sk), PublicKey(pk)))
        Ok((sk, pk))
    }

    pub fn encode<const P: u64>(m: &Rq<P, 1>) -> T64 {
        let delta = u64::MAX / P; // floored
        let coeffs = m.coeffs();
        // Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
        T64(coeffs[0].0 * delta)
    }
    pub fn decode<const P: u64>(p: &T64) -> Rq<P, 1> {
        let p = p.mul_div_round(P, u64::MAX);
        Rq::<P, 1>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(rng: impl Rng, sk: &SecretKey<K>, p: &T64) -> Result<Self> {
        let glwe = GLWE::encrypt_s(rng, &sk, p)?;
        Ok(Self(glwe))
    }
    pub fn encrypt(rng: impl Rng, pk: &PublicKey<K>, p: &T64) -> Result<Self> {
        let glwe = GLWE::encrypt(rng, &pk, p)?;
        Ok(Self(glwe))
    }
    pub fn decrypt(&self, sk: &SecretKey<K>) -> T64 {
        self.0.decrypt(&sk)
    }
}

impl<const K: usize> Add<TLWE<K>> for TLWE<K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}
impl<const K: usize> AddAssign for TLWE<K> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}
impl<const K: usize> Sum<TLWE<K>> for TLWE<K> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = TLWE::<K>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<const K: usize> Sub<TLWE<K>> for TLWE<K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

// plaintext addition
impl<const K: usize> Add<T64> for TLWE<K> {
    type Output = Self;
    fn add(self, plaintext: T64) -> Self {
        let a: TR<T64, K> = self.0 .0;
        let b: T64 = self.0 .1 + plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext substraction
impl<const K: usize> Sub<T64> for TLWE<K> {
    type Output = Self;
    fn sub(self, plaintext: T64) -> Self {
        let a: TR<T64, K> = self.0 .0;
        let b: T64 = self.0 .1 - plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext multiplication
impl<const K: usize> Mul<T64> for TLWE<K> {
    type Output = Self;
    fn mul(self, plaintext: T64) -> Self {
        let a: TR<T64, K> = TR(self.0 .0 .0.iter().map(|r_i| *r_i * plaintext).collect());
        let b: T64 = self.0 .1 * plaintext;
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
        const K: usize = 16;
        type S = TLWE<K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            dbg!(&m);
            let p: T64 = S::encode::<T>(&m);
            dbg!(&p);

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
        const K: usize = 16;
        type S = TLWE<K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let p1: T64 = S::encode::<T>(&m1); // plaintext
            let p2: T64 = S::encode::<T>(&m2); // plaintext

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
        const K: usize = 16;
        type S = TLWE<K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let p1: T64 = S::encode::<T>(&m1); // plaintext
            let p2: T64 = S::encode::<T>(&m2); // plaintext

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
        const K: usize = 16;
        type S = TLWE<K>;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let m1 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let p1: T64 = S::encode::<T>(&m1);
            // don't scale up p2, set it directly from m2
            // let p2: T64 = Tn(array::from_fn(|i| T64(m2.coeffs()[i].0)));
            let p2: T64 = T64(m2.coeffs()[0].0);

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: T64 = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), m3_recovered);
        }

        Ok(())
    }
}
