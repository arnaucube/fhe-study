use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, RingParam, Rq, Tn, T64, TR};

use crate::tlev::TLev;
use crate::{
    tglwe::TGLWE,
    tlwe::{PublicKey, SecretKey, TLWE},
};
use gfhe::glwe::{Param, GLWE};

/// vector of length K+1 = [K], [1]
#[derive(Clone, Debug)]
pub struct TGSW(pub(crate) Vec<TLev>, TLev);

impl TGSW {
    pub fn encrypt_s(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        sk: &SecretKey,
        m: &T64,
    ) -> Result<Self> {
        let a: Vec<TLev> = (0..param.k)
            .map(|i| TLev::encrypt_s(&mut rng, &param, beta, l, sk, &(-sk.0 .0.r[i] * *m)))
            .collect::<Result<Vec<_>>>()?;
        let b: TLev = TLev::encrypt_s(&mut rng, &param, beta, l, sk, m)?;
        Ok(Self(a, b))
    }

    pub fn decrypt(&self, sk: &SecretKey, beta: u32) -> T64 {
        self.1.decrypt(sk, beta)
    }
    pub fn from_tlwe(_tlwe: TLWE) -> Self {
        todo!()
    }

    pub fn cmux(bit: Self, ct1: TLWE, ct2: TLWE) -> TLWE {
        ct1.clone() + (bit * (ct2 - ct1))
    }
}

/// External product TGSW x TLWE
impl Mul<TLWE> for TGSW {
    type Output = TLWE;

    fn mul(self, tlwe: TLWE) -> TLWE {
        let beta: u32 = 2;
        let l: u32 = 64; // TODO wip

        // since N=1, each tlwe element is a vector of length=1, decomposed into
        // l elements, and we have K of them
        let tlwe_ab: Vec<T64> = [tlwe.0 .0.r.clone(), vec![tlwe.0 .1]].concat();

        let tgsw_ab: Vec<TLev> = [self.0.clone(), vec![self.1]].concat();
        assert_eq!(tgsw_ab.len(), tlwe_ab.len());

        let r: TLWE = zip_eq(tgsw_ab, tlwe_ab)
            .map(|(tlev_i, tlwe_i)| tlev_i * tlwe_i.decompose(beta, l))
            .sum();
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
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 2, // plaintext modulus
        };
        let beta: u32 = 2;
        let l: u32 = 16;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..50 {
            let (sk, _) = TLWE::new_key(&mut rng, &param)?;

            let m: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p: T64 = TLev::encode(&param, &m); // plaintext

            let c = TGSW::encrypt_s(&mut rng, &param, beta, l, &sk, &p)?;
            let p_recovered = c.decrypt(&sk, beta);
            let m_recovered = TLev::decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_external_product() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 32,
            t: 2, // plaintext modulus
        };
        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..50 {
            let (sk, _) = TLWE::new_key(&mut rng, &param)?;

            let m1: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLev::encode(&param, &m1);

            let m2: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p2: T64 = TLWE::encode(&param, &m2); // scaled by delta

            let tgsw = TGSW::encrypt_s(&mut rng, &param, beta, l, &sk, &p1)?;
            let tlwe = TLWE::encrypt_s(&mut rng, &param, &sk, &p2)?;

            let res: TLWE = tgsw * tlwe;

            // let p_recovered = res.decrypt(&sk, beta);
            let p_recovered = res.decrypt(&sk);
            // downscaled by delta^-1
            let res_recovered = TLWE::decode(&param, &p_recovered);

            // assert_eq!(m1 * m2, m_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), res_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_cmux() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 32,
            t: 2, // plaintext modulus
        };

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..50 {
            let (sk, _) = TLWE::new_key(&mut rng, &param)?;

            let m1: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLWE::encode(&param, &m1); // scaled by delta

            let m2: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p2: T64 = TLWE::encode(&param, &m2); // scaled by delta

            for bit_raw in 0..2 {
                let bit = TGSW::encrypt_s(&mut rng, &param, beta, l, &sk, &T64(bit_raw))?;

                let c1 = TLWE::encrypt_s(&mut rng, &param, &sk, &p1)?;
                let c2 = TLWE::encrypt_s(&mut rng, &param, &sk, &p2)?;

                let res: TLWE = TGSW::cmux(bit, c1, c2);

                let p_recovered = res.decrypt(&sk);
                // downscaled by delta^-1
                let res_recovered = TLWE::decode(&param, &p_recovered);

                if bit_raw == 0 {
                    assert_eq!(m1, res_recovered);
                } else {
                    assert_eq!(m2, res_recovered);
                }
            }
        }

        Ok(())
    }
}
