//! Implementation of TFHE https://eprint.iacr.org/2018/421.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::array;

use arith::{Ring, Rq, Tn, T64};
use gfhe::{glwe, GLWE};

pub mod tlev;
pub mod tlwe;

#[derive(Clone, Debug)]
pub struct SecretKey<const K: usize>(glwe::SecretKey<Tn<1>, K>);
#[derive(Clone, Debug)]
pub struct PublicKey<const K: usize>(glwe::PublicKey<Tn<1>, K>);

#[derive(Clone, Debug)]
pub struct TLWE<const K: usize>(pub GLWE<Tn<1>, K>);

impl<const K: usize> TLWE<K> {
    pub fn new_key(rng: impl Rng) -> Result<(SecretKey<K>, PublicKey<K>)> {
        let (sk, pk) = GLWE::new_key(rng)?;
        Ok((SecretKey(sk), PublicKey(pk)))
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
    pub fn encrypt_s(rng: impl Rng, sk: &SecretKey<K>, p: &Tn<1>) -> Result<Self> {
        let glwe = GLWE::encrypt_s(rng, &sk.0, p)?;
        Ok(Self(glwe))
    }

    pub fn encrypt(rng: impl Rng, pk: &PublicKey<K>, p: &Tn<1>) -> Result<Self> {
        let glwe = GLWE::encrypt(rng, &pk.0, p)?;
        Ok(Self(glwe))
    }

    pub fn decrypt(&self, sk: &SecretKey<K>) -> Tn<1> {
        self.0.decrypt(&sk.0)
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        const T: u64 = 128; // plaintext modulus
        const K: usize = 16;
        type S = TLWE<K>;

        // let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_f64, T as f64);
            let m = Rq::<T, 1>::rand(&mut rng, msg_dist); // msg

            // let m: Rq<Q, N> = m.remodule::<Q>();

            let p = S::encode::<T>(&m); // plaintext
            let c = S::encrypt(&mut rng, &pk, &p)?; // ciphertext
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }
}
