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

use arithmetic::PR;

// error deviation for the Gaussian(Normal) distribution
// sigma=3.2 from: https://eprint.iacr.org/2022/162.pdf page 5
const ERR_SIGMA: f64 = 3.2;

#[derive(Clone, Debug)]
pub struct SecretKey<const Q: u64, const N: usize>(PR<Q, N>);

#[derive(Clone, Debug)]
pub struct PublicKey<const Q: u64, const N: usize>(PR<Q, N>, PR<Q, N>);

// RLWE ciphertext
#[derive(Clone, Debug)]
pub struct RLWE<const Q: u64, const N: usize>(PR<Q, N>, PR<Q, N>);

impl<const Q: u64, const N: usize> RLWE<Q, N> {
    fn add(lhs: Self, rhs: Self) -> Self {
        RLWE::<Q, N>(lhs.0 + rhs.0, lhs.1 + rhs.1)
    }
    fn mul(lhs: Self, rhs: Self) -> Self {
        todo!()
    }
}

impl<const Q: u64, const N: usize> ops::Add<RLWE<Q, N>> for RLWE<Q, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::add(self, rhs)
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
        let s = PR::<Q, N>::rand_u64(&mut rng, Xi_key)?;

        // pk = (-a * s + e, a)
        let a = PR::<Q, N>::rand_u64(&mut rng, Uniform::new(0_u64, Q))?;
        let e = PR::<Q, N>::rand_f64(&mut rng, Xi_err)?;
        let pk: PublicKey<Q, N> = PublicKey((&(-a) * &s) + e, a.clone());
        Ok((SecretKey(s), pk))
    }

    pub fn encrypt(mut rng: impl Rng, pk: &PublicKey<Q, N>, m: &PR<T, N>) -> Result<RLWE<Q, N>> {
        let Xi_key = Uniform::new(-1_f64, 1_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u = PR::<Q, N>::rand_f64(&mut rng, Xi_key)?;
        let e_1 = PR::<Q, N>::rand_f64(&mut rng, Xi_err)?;
        let e_2 = PR::<Q, N>::rand_f64(&mut rng, Xi_err)?;

        // migrate m's coeffs to the bigger modulus Q (from T)
        let m = PR::<Q, N>::from_vec_u64(m.coeffs().iter().map(|m_i| m_i.0).collect());
        let c0 = &pk.0 * &u + e_1 + m * Self::DELTA;
        let c1 = &pk.1 * &u + e_2;
        Ok(RLWE::<Q, N>(c0, c1))
    }

    pub fn decrypt(sk: &SecretKey<Q, N>, c: &RLWE<Q, N>) -> PR<T, N> {
        let cs = c.0 + c.1 * sk.0; // done in mod q
        let r: Vec<u64> = cs
            .coeffs()
            .iter()
            .map(|e| ((T as f64 * e.0 as f64) / Q as f64).round() as u64)
            .collect();
        PR::<T, N>::from_vec_u64(r)
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

        let (sk, pk) = S::new_key(&mut rng)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m = PR::<T, N>::rand_u64(&mut rng, msg_dist)?;

        let c = S::encrypt(rng, &pk, &m)?;
        let m_recovered = S::decrypt(&sk, &c);

        assert_eq!(m, m_recovered);

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 32;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m1 = PR::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let m2 = PR::<T, N>::rand_u64(&mut rng, msg_dist)?;

        let c1 = S::encrypt(&mut rng, &pk, &m1)?;
        let c2 = S::encrypt(&mut rng, &pk, &m2)?;

        let c3 = c1 + c2;

        let m3_recovered = S::decrypt(&sk, &c3);

        assert_eq!(m1 + m2, m3_recovered);

        Ok(())
    }
}
