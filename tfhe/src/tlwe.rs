use anyhow::Result;
use itertools::zip_eq;
use rand::distributions::Standard;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::array;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, Rq, Tn, T64, TR};

const ERR_SIGMA: f64 = 3.2;

#[derive(Clone, Debug)]
pub struct TLWE<const K: usize>(TR<Tn<1>, K>, Tn<1>);

#[derive(Clone, Debug)]
pub struct SecretKey<const K: usize>(TR<Tn<1>, K>);
#[derive(Clone, Debug)]
pub struct PublicKey<const K: usize>(Tn<1>, TR<Tn<1>, K>);

impl<const K: usize> TLWE<K> {
    pub fn zero() -> Self {
        Self(TR::zero(), Tn::zero())
    }

    pub fn new_key(mut rng: impl Rng) -> Result<(SecretKey<K>, PublicKey<K>)> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let s: TR<Tn<1>, K> = TR::rand(&mut rng, Xi_key);
        let a: TR<Tn<1>, K> = TR::rand(&mut rng, Standard);
        let e = Tn::rand(&mut rng, Xi_err);

        let pk: PublicKey<K> = PublicKey((&a * &s) + e, a);
        Ok((SecretKey(s), pk))
    }

    pub fn encode<const P: u64>(m: &Rq<P, 1>) -> Tn<1> {
        let delta = u64::MAX / P; // floored
        let coeffs = m.coeffs();
        Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
    }
    pub fn decode<const P: u64>(p: &Tn<1>) -> Rq<P, 1> {
        let p = p.mul_div_round(P, u64::MAX);
        Rq::<P, 1>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(mut rng: impl Rng, sk: &SecretKey<K>, m: &Tn<1>) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<Tn<1>, K> = TR::rand(&mut rng, Xi_key);
        let e = Tn::rand(&mut rng, Xi_err);

        let b: Tn<1> = (&a * &sk.0) + *m + e;
        Ok(Self(a, b))
    }
    pub fn encrypt(mut rng: impl Rng, pk: &PublicKey<K>, m: &Tn<1>) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u: Tn<1> = Tn::rand(&mut rng, Xi_key);

        let e0: Tn<1> = Tn::rand(&mut rng, Xi_err);
        let e1 = TR::<Tn<1>, K>::rand(&mut rng, Xi_err);

        let b: Tn<1> = pk.0 * u + *m + e0;
        let d: TR<Tn<1>, K> = &pk.1 * &u + e1;

        Ok(Self(d, b))
    }
    pub fn decrypt(&self, sk: &SecretKey<K>) -> Tn<1> {
        let (d, b): (TR<Tn<1>, K>, Tn<1>) = (self.0.clone(), self.1);
        b - &d * &sk.0
    }
}

impl<const K: usize> Add<TLWE<K>> for TLWE<K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let a: TR<Tn<1>, K> = self.0 + other.0;
        let b: Tn<1> = self.1 + other.1;
        Self(a, b)
    }
}
impl<const K: usize> AddAssign for TLWE<K> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..K {
            self.0 .0[i] = self.0 .0[i] + rhs.0 .0[i];
        }
        self.1 = self.1 + rhs.1;
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
        let a: TR<Tn<1>, K> = self.0 - other.0;
        let b: Tn<1> = self.1 - other.1;
        Self(a, b)
    }
}

// plaintext addition
impl<const K: usize> Add<Tn<1>> for TLWE<K> {
    type Output = Self;
    fn add(self, plaintext: Tn<1>) -> Self {
        let a: TR<Tn<1>, K> = self.0;
        let b: Tn<1> = self.1 + plaintext;
        Self(a, b)
    }
}
// plaintext substraction
impl<const K: usize> Sub<Tn<1>> for TLWE<K> {
    type Output = Self;
    fn sub(self, plaintext: Tn<1>) -> Self {
        let a: TR<Tn<1>, K> = self.0;
        let b: Tn<1> = self.1 - plaintext;
        Self(a, b)
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

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            dbg!(&m);
            let p: Tn<1> = S::encode::<T>(&m);
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

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m1 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<1> = S::encode::<T>(&m1); // plaintext
            let p2: Tn<1> = S::encode::<T>(&m2); // plaintext

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

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m1 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, 1>::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<1> = S::encode::<T>(&m1); // plaintext
            let p2: Tn<1> = S::encode::<T>(&m2); // plaintext

            let c1 = S::encrypt(&mut rng, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = S::decode::<T>(&p3_recovered);

            assert_eq!((m1 + m2).remodule::<T>(), m3_recovered);
        }

        Ok(())
    }
}
