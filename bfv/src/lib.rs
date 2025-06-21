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

use arithmetic::{Zq, PR};

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
}

impl<const Q: u64, const N: usize> ops::Add<RLWE<Q, N>> for RLWE<Q, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::add(self, rhs)
    }
}

impl<const Q: u64, const N: usize, const T: u64> ops::Add<&PR<T, N>> for &RLWE<Q, N> {
    type Output = RLWE<Q, N>;
    fn add(self, rhs: &PR<T, N>) -> Self::Output {
        // todo!()
        BFV::<Q, N, T>::add_const(self, rhs)
    }
}
impl<const Q: u64, const N: usize, const T: u64> ops::Mul<&PR<T, N>> for &RLWE<Q, N> {
    type Output = RLWE<Q, N>;
    fn mul(self, rhs: &PR<T, N>) -> Self::Output {
        // todo!()
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

    fn add_const(c: &RLWE<Q, N>, m: &PR<T, N>) -> RLWE<Q, N> {
        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule::<Q>();
        RLWE::<Q, N>(c.0 + m * Self::DELTA, c.1)
    }
    fn mul_const(c: &RLWE<Q, N>, m: &PR<T, N>) -> RLWE<Q, N> {
        // assuming T<Q, move m from Zq<T> to Zq<Q>
        let m = m.remodule::<Q>();
        RLWE::<Q, N>(c.0 * m * Self::DELTA, c.1)
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

    #[test]
    fn test_constant_add_mul() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 32;
        const T: u64 = 4; // plaintext modulus
        type S = BFV<Q, N, T>;

        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m1 = PR::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let m2_const = PR::<T, N>::rand_u64(&mut rng, msg_dist)?;

        let c1 = S::encrypt(&mut rng, &pk, &m1)?;

        let c3_add = &c1 + &m2_const;
        let c3_mul = &c1 * &m2_const;

        let m3_add_recovered = S::decrypt(&sk, &c3_add);
        let m3_mul_recovered = S::decrypt(&sk, &c3_mul);

        assert_eq!(m1 + m2_const, m3_add_recovered);

        let mut mul_res = naive_poly_mul::<T>(&m1.coeffs().to_vec(), &m2_const.coeffs().to_vec());
        arithmetic::ring::modulus::<T, N>(&mut mul_res);
        dbg!(&mul_res);
        let mul_res_2 =
            naive_poly_mul_2::<T, N>(&m1.coeffs().to_vec(), &m2_const.coeffs().to_vec());
        assert_eq!(mul_res, mul_res_2);
        let mul_res = PR::<T, N>::from_vec(mul_res);
        assert_eq!(mul_res.coeffs(), m3_mul_recovered.coeffs());

        Ok(())
    }

    fn naive_poly_mul<const T: u64>(a: &[Zq<T>], b: &[Zq<T>]) -> Vec<Zq<T>> {
        let mut result: Vec<Zq<T>> = vec![Zq::zero(); a.len() + b.len() - 1];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                result[i + j] = result[i + j] + (ai * bj);
            }
        }
        result
    }
    fn naive_poly_mul_2<const T: u64, const N: usize>(
        poly1: &[Zq<T>],
        poly2: &[Zq<T>],
    ) -> Vec<Zq<T>> {
        let degree1 = poly1.len();
        let degree2 = poly2.len();

        // The degree of the resulting polynomial will be degree1 + degree2 - 1
        let mut result = vec![Zq::zero(); degree1 + degree2 - 1];

        // Perform the multiplication
        for i in 0..degree1 {
            for j in 0..degree2 {
                result[i + j] = result[i + j] + poly1[i] * poly2[j];
            }
        }

        // Reduce the result modulo x^N + 1
        let mut reduced_result = vec![Zq::zero(); N];

        for i in 0..result.len() {
            let mod_index = i % N; // wrap around using modulo N
            reduced_result[mod_index] += result[i];
        }

        // Return the reduced polynomial
        reduced_result
    }
}
