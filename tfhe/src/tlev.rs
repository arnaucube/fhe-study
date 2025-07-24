use anyhow::Result;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, Rq, Tn, T64, TR};

use crate::tlwe::{PublicKey, SecretKey, TLWE};

#[derive(Clone, Debug)]
pub struct TLev<const K: usize>(pub(crate) Vec<TLWE<K>>);

impl<const K: usize> TLev<K> {
    pub fn encode<const T: u64>(m: &Rq<T, 1>) -> Tn<1> {
        let coeffs = m.coeffs();
        Tn(array::from_fn(|i| T64(coeffs[i].0)))
    }
    pub fn decode<const T: u64>(p: &Tn<1>) -> Rq<T, 1> {
        Rq::<T, 1>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<K>,
        m: &Tn<1>,
    ) -> Result<Self> {
        let tlev: Vec<TLWE<K>> = (1..l + 1)
            .map(|i| {
                TLWE::<K>::encrypt(&mut rng, pk, &(*m * (u64::MAX / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<K>,
        m: &Tn<1>,
    ) -> Result<Self> {
        let tlev: Vec<TLWE<K>> = (1..l + 1)
            .map(|i| {
                TLWE::<K>::encrypt_s(&mut rng, sk, &(*m * (u64::MAX / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }

    pub fn decrypt(&self, sk: &SecretKey<K>, beta: u32) -> Tn<1> {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, u64::MAX)
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        const T: u64 = 2; // plaintext modulus
        const K: usize = 16;
        type S = TLev<K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = TLWE::<K>::new_key(&mut rng)?;

            let m: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p: Tn<1> = S::encode::<T>(&m); // plaintext

            let c = S::encrypt(&mut rng, beta, l, &pk, &p)?;
            let p_recovered = c.decrypt(&sk, beta);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }
}
