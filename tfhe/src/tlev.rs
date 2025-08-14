use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::ops::Mul;

use arith::{Ring, RingParam, Rq, T64};

use crate::tlwe::{PublicKey, SecretKey, TLWE};
use gfhe::glwe::Param;

#[derive(Clone, Debug)]
pub struct TLev(pub(crate) Vec<TLWE>);

impl TLev {
    pub fn encode(param: &Param, m: &Rq) -> T64 {
        assert_eq!(m.param.n, 1);
        assert_eq!(param.t, m.param.q);

        let coeffs = m.coeffs();
        T64(coeffs[0].v) // N=1, so take the only coeff
    }
    pub fn decode(param: &Param, p: &T64) -> Rq {
        Rq::from_vec_u64(
            &RingParam { q: param.t, n: 1 },
            p.coeffs().iter().map(|c| c.0).collect(),
        )
    }
    pub fn encrypt(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        pk: &PublicKey,
        m: &T64,
    ) -> Result<Self> {
        debug_assert_eq!(pk.1.k, param.k);

        let tlev: Vec<TLWE> = (1..l as u64 + 1)
            .map(|i| {
                let aux = if i < 64 {
                    *m * (u64::MAX / (1u64 << i))
                } else {
                    // 1<<64 would overflow, and anyways we're dividing u64::MAX
                    // by it, which would be equal to 1
                    *m
                };
                TLWE::encrypt(&mut rng, param, pk, &aux)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        param: &Param,
        _beta: u32, // TODO rm, and make beta=2 always
        l: u32,
        sk: &SecretKey,
        m: &T64,
    ) -> Result<Self> {
        debug_assert_eq!(sk.0 .0.k, param.k);

        let tlev: Vec<TLWE> = (1..l as u64 + 1)
            .map(|i| {
                let aux = if i < 64 {
                    *m * (u64::MAX / (1u64 << i))
                } else {
                    // 1<<64 would overflow, and anyways we're dividing u64::MAX
                    // by it, which would be equal to 1
                    *m
                };
                TLWE::encrypt_s(&mut rng, &param, sk, &aux)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }

    pub fn decrypt(&self, sk: &SecretKey, beta: u32) -> T64 {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, u64::MAX)
    }
}
// TODO review u64::MAX, since is -1 of the value we actually want

impl TLev {
    pub fn iter(&self) -> std::slice::Iter<TLWE> {
        self.0.iter()
    }
}

// dot product between a TLev and Vec<T64>, usually Vec<T64> comes from a
// decomposition of T64
// TLev * Vec<T64> --> TLWE
impl Mul<Vec<T64>> for TLev {
    type Output = TLWE;
    fn mul(self, v: Vec<T64>) -> Self::Output {
        assert_eq!(self.0.len(), v.len());

        // l TLWES
        let tlwes: Vec<TLWE> = self.0;
        let r: TLWE = zip_eq(v, tlwes).map(|(a_d_i, glwe_i)| glwe_i * a_d_i).sum();
        r
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 2, // plaintext modulus
        };

        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p: T64 = TLev::encode(&param, &m); // plaintext

            let c = TLev::encrypt(&mut rng, &param, beta, l, &pk, &p)?;
            let p_recovered = c.decrypt(&sk, beta);
            let m_recovered = TLev::decode(&param, &p_recovered);

            assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_tlev_vect64_product() -> Result<()> {
        let param = Param {
            err_sigma: 0.1, // WIP
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 2, // plaintext modulus
        };

        let beta: u32 = 2;
        // let l: u32 = 16;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m1: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLev::encode(&param, &m1);
            let p2: T64 = TLev::encode(&param, &m2);

            let c1 = TLev::encrypt(&mut rng, &param, beta, l, &pk, &p1)?;
            let c2 = p2.decompose(beta, l);

            let c3 = c1 * c2;

            let p_recovered = c3.decrypt(&sk);
            let m_recovered = TLev::decode(&param, &p_recovered);

            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), m_recovered);
        }

        Ok(())
    }
}
