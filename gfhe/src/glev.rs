use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::ops::{Add, Mul};

use arith::{Ring, TR};

use crate::glwe::{PublicKey, SecretKey, GLWE};

const ERR_SIGMA: f64 = 3.2;

// l GLWEs
#[derive(Clone, Debug)]
pub struct GLev<R: Ring, const K: usize>(pub(crate) Vec<GLWE<R, K>>);

impl<R: Ring, const K: usize> GLev<R, K> {
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<R, K>,
        m: &R,
    ) -> Result<Self> {
        let glev: Vec<GLWE<R, K>> = (0..l)
            .map(|i| {
                GLWE::<R, K>::encrypt(&mut rng, pk, &(*m * (R::Q / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<R, K>,
        m: &R,
        // delta: u64,
    ) -> Result<Self> {
        let glev: Vec<GLWE<R, K>> = (1..l + 1)
            .map(|i| {
                GLWE::<R, K>::encrypt_s(&mut rng, sk, &(*m * (R::Q / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }

    pub fn decrypt<const T: u64>(&self, sk: &SecretKey<R, K>, beta: u32) -> R {
        let pt = self.0[1].decrypt(sk);
        pt.mul_div_round(beta as u64, R::Q)
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;
    use arith::Rq;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 2; // plaintext modulus
        const K: usize = 16;
        type S = GLev<Rq<Q, N>, K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        // let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..200 {
            let (sk, pk) = GLWE::<Rq<Q, N>, K>::new_key(&mut rng)?;

            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m: Rq<Q, N> = m.remodule::<Q>();

            let c = S::encrypt(&mut rng, beta, l, &pk, &m)?;
            let m_recovered = c.decrypt::<T>(&sk, beta);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }
}
