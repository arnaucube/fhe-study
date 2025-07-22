use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::ops::{Add, Mul};

use arith::{Ring, Rq, TR};

use crate::glwe::{PublicKey, SecretKey, GLWE};

const ERR_SIGMA: f64 = 3.2;

// l GLWEs
#[derive(Clone, Debug)]
pub struct GLev<const Q: u64, const N: usize, const K: usize>(pub(crate) Vec<GLWE<Q, N, K>>);

impl<const Q: u64, const N: usize, const K: usize> GLev<Q, N, K> {
    pub fn encode<const T: u64>(m: &Rq<T, N>) -> Rq<Q, N> {
        m.remodule::<Q>()
    }
    pub fn decode<const T: u64>(p: &Rq<Q, N>) -> Rq<T, N> {
        p.remodule::<T>()
    }
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<Q, N, K>,
        m: &Rq<Q, N>,
    ) -> Result<Self> {
        let glev: Vec<GLWE<Q, N, K>> = (1..l + 1)
            .map(|i| {
                GLWE::<Q, N, K>::encrypt(&mut rng, pk, &(*m * (Q / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<Q, N, K>,
        m: &Rq<Q, N>,
    ) -> Result<Self> {
        let glev: Vec<GLWE<Q, N, K>> = (1..l + 1)
            .map(|i| {
                GLWE::<Q, N, K>::encrypt_s(&mut rng, sk, &(*m * (Q / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }

    pub fn decrypt<const T: u64>(&self, sk: &SecretKey<Q, N, K>, beta: u32) -> Rq<Q, N> {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, Q)
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
        const N: usize = 128;
        const T: u64 = 2; // plaintext modulus
        const K: usize = 16;
        type S = GLev<Q, N, K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = GLWE::<Q, N, K>::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let p: Rq<Q, N> = S::encode::<T>(&m); // plaintext

            let c = S::encrypt(&mut rng, beta, l, &pk, &p)?;
            let p_recovered = c.decrypt::<T>(&sk, beta);
            let m_recovered = S::decode::<T>(&p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }
}
