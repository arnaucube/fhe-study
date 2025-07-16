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
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<Q, N, K>,
        m: &Rq<Q, N>,
        // delta: u64,
    ) -> Result<Self> {
        let glev: Vec<GLWE<Q, N, K>> = (0..l)
            .map(|i| {
                GLWE::<Q, N, K>::encrypt(&mut rng, pk, &(*m * (Q / beta.pow(i as u32) as u64)), 1)
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
        // delta: u64,
    ) -> Result<Self> {
        let glev: Vec<GLWE<Q, N, K>> = (0..l)
            .map(|i| {
                GLWE::<Q, N, K>::encrypt_s(&mut rng, sk, &(*m * (Q / beta.pow(i as u32) as u64)), 1)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }

    pub fn decrypt<const T: u64>(&self, sk: &SecretKey<Q, N, K>, delta: u64) -> Rq<Q, N> {
        self.0[1].decrypt::<T>(sk, delta)
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

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = GLWE::<Q, N, K>::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m: Rq<Q, N> = m.remodule::<Q>();

            let c = S::encrypt(&mut rng, beta, l, &pk, &m)?;
            let m_recovered = c.decrypt::<T>(&sk, delta);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }
}
