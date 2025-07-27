use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, Rq, Tn, T64, TR};

use crate::tlev::TLev;
use crate::tlwe::{PublicKey, SecretKey, TLWE};
use gfhe::glwe::GLWE;

/// vector of length K+1 = [K], [1]
#[derive(Clone, Debug)]
pub struct TGSW<const K: usize>(pub(crate) Vec<TLev<K>>, TLev<K>);

impl<const K: usize> TGSW<K> {
    pub fn encrypt_s(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<K>,
        m: &T64,
    ) -> Result<Self> {
        let a: Vec<TLev<K>> = (0..K)
            .map(|i| TLev::encrypt_s(&mut rng, beta, l, sk, &(-sk.0 .0[i] * *m)))
            .collect::<Result<Vec<_>>>()?;
        let b: TLev<K> = TLev::encrypt_s(&mut rng, beta, l, sk, m)?;
        Ok(Self(a, b))
    }

    pub fn decrypt(&self, sk: &SecretKey<K>, beta: u32) -> T64 {
        self.1.decrypt(sk, beta)
    }
    pub fn from_tlwe(_tlwe: TLWE<K>) -> Self {
        todo!()
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
        type S = TGSW<K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..50 {
            let (sk, _) = TLWE::<K>::new_key(&mut rng)?;

            let m: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p: T64 = TLev::<K>::encode::<T>(&m); // plaintext

            let c = S::encrypt_s(&mut rng, beta, l, &sk, &p)?;
            let p_recovered = c.decrypt(&sk, beta);
            let m_recovered = TLev::<K>::decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }
}
