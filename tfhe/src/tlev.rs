use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, Rq, Tn, T64, TR};

use crate::tlwe::{PublicKey, SecretKey, TLWE};

#[derive(Clone, Debug)]
pub struct TLev<const K: usize>(pub(crate) Vec<TLWE<K>>);

impl<const K: usize> TLev<K> {
    pub fn encode<const T: u64>(m: &Rq<T, 1>) -> T64 {
        let coeffs = m.coeffs();
        T64(coeffs[0].0) // N=1, so take the only coeff
    }
    pub fn decode<const T: u64>(p: &T64) -> Rq<T, 1> {
        Rq::<T, 1>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<K>,
        m: &T64,
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
        _beta: u32, // TODO rm, and make beta=2 always
        l: u32,
        sk: &SecretKey<K>,
        m: &T64,
    ) -> Result<Self> {
        let tlev: Vec<TLWE<K>> = (1..l as u64 + 1)
            .map(|i| {
                let aux = if i < 64 {
                    *m * (u64::MAX / (1u64 << i))
                } else {
                    // 1<<64 would overflow, and anyways we're dividing u64::MAX
                    // by it, which would be equal to 1
                    *m
                };
                TLWE::<K>::encrypt_s(&mut rng, sk, &aux)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }

    pub fn decrypt(&self, sk: &SecretKey<K>, beta: u32) -> T64 {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, u64::MAX)
    }
}
// TODO review u64::MAX, since is -1 of the value we actually want

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
            let p: T64 = S::encode::<T>(&m); // plaintext

            let c = S::encrypt(&mut rng, beta, l, &pk, &p)?;
            let p_recovered = c.decrypt(&sk, beta);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }
}
