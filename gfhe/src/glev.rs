use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::ops::Mul;

use arith::Ring;

use crate::glwe::{Param, PublicKey, SecretKey, GLWE};

// l GLWEs
#[derive(Clone, Debug)]
pub struct GLev<R: Ring>(pub(crate) Vec<GLWE<R>>);

impl<R: Ring> GLev<R> {
    pub fn encrypt(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        pk: &PublicKey<R>,
        m: &R,
    ) -> Result<Self> {
        let glev: Vec<GLWE<R>> = (0..l)
            .map(|i| {
                GLWE::<R>::encrypt(
                    &mut rng,
                    param,
                    pk,
                    &(m.clone() * (param.ring.q / beta.pow(i as u32) as u64)),
                )
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        sk: &SecretKey<R>,
        m: &R,
    ) -> Result<Self> {
        let glev: Vec<GLWE<R>> = (1..l + 1)
            .map(|i| {
                GLWE::<R>::encrypt_s(
                    &mut rng,
                    param,
                    sk,
                    &(m.clone() * (param.ring.q / beta.pow(i as u32) as u64)), // TODO rm clone
                )
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(glev))
    }

    pub fn decrypt(&self, param: &Param, sk: &SecretKey<R>, beta: u32) -> R {
        let pt = self.0[1].decrypt(sk);
        pt.mul_div_round(beta as u64, param.ring.q)
    }
}

// dot product between a GLev and Vec<R>.
// Used for operating decompositions with KSK_i.
// GLev * Vec<R> --> GLWE
impl<R: Ring> Mul<Vec<R>> for GLev<R> {
    type Output = GLWE<R>;
    fn mul(self, v: Vec<R>) -> GLWE<R> {
        debug_assert_eq!(self.0.len(), v.len());
        // TODO debug_assert_eq of param

        // l times GLWES
        let glwes: Vec<GLWE<R>> = self.0;

        // l iterations
        let r: GLWE<R> = zip_eq(v, glwes).map(|(v_i, glwe_i)| glwe_i * v_i).sum();
        r
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;
    use arith::{RingParam, Rq};

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        let param = Param {
            err_sigma: crate::glwe::ERR_SIGMA,
            ring: RingParam {
                q: 2u64.pow(16) + 1,
                n: 128,
            },
            k: 16,
            t: 2, // plaintext modulus
        };
        type S = GLev<Rq>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = GLWE::<Rq>::new_key(&mut rng, &param)?;

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m: Rq = m.remodule(param.ring.q);

            let c = S::encrypt(&mut rng, &param, beta, l, &pk, &m)?;
            let m_recovered = c.decrypt(&param, &sk, beta);

            assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));
        }

        Ok(())
    }
}
