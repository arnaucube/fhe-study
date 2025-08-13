use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, RingParam, Rq, Tn, T64, TR};

use crate::tglwe::{PublicKey, SecretKey, TGLWE};
use gfhe::glwe::{Param, GLWE};

/// vector of length K+1 = ([K * TGLev], [1 * TGLev])
#[derive(Clone, Debug)]
pub struct TGGSW(pub(crate) Vec<TGLev>, TGLev);

impl TGGSW {
    pub fn encrypt_s(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        sk: &SecretKey,
        m: &Tn,
    ) -> Result<Self> {
        debug_assert_eq!(sk.0 .0.k, param.k);

        let a: Vec<TGLev> = (0..param.k)
            .map(|i| TGLev::encrypt_s(&mut rng, param, beta, l, sk, &(&-sk.0 .0.r[i].clone() * m)))
            // TODO rm clone
            .collect::<Result<Vec<_>>>()?;
        let b: TGLev = TGLev::encrypt_s(&mut rng, &param, beta, l, sk, m)?;
        Ok(Self(a, b))
    }

    pub fn decrypt(&self, sk: &SecretKey, beta: u32) -> Tn {
        self.1.decrypt(sk, beta)
    }

    pub fn cmux(bit: Self, ct1: TGLWE, ct2: TGLWE) -> TGLWE {
        ct1.clone() + (bit * (ct2 - ct1))
    }
}

/// External product tggsw x tglwe
impl Mul<TGLWE> for TGGSW {
    type Output = TGLWE;

    fn mul(self, tglwe: TGLWE) -> TGLWE {
        let beta: u32 = 2;
        let l: u32 = 64; // TODO wip

        let tglwe_ab: Vec<Tn> = [tglwe.0 .0.r.clone(), vec![tglwe.0 .1]].concat();

        let tgsw_ab: Vec<TGLev> = [self.0.clone(), vec![self.1]].concat();
        assert_eq!(tgsw_ab.len(), tglwe_ab.len());

        let r: TGLWE = zip_eq(tgsw_ab, tglwe_ab)
            .map(|(tlev_i, tglwe_i)| tlev_i * tglwe_i.decompose(beta, l))
            .sum();
        r
    }
}

#[derive(Clone, Debug)]
pub struct TGLev(pub(crate) Vec<TGLWE>);

impl TGLev {
    pub fn encode(param: &Param, m: &Rq) -> Tn {
        debug_assert_eq!(param.t, m.param.q); // plaintext modulus

        Tn {
            param: param.ring,
            coeffs: m.coeffs().iter().map(|c_i| T64(c_i.v)).collect(),
        }
    }
    pub fn decode(param: &Param, p: &Tn) -> Rq {
        Rq::from_vec_u64(&param.pt(), p.coeffs().iter().map(|c| c.0).collect())
    }
    pub fn encrypt(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        pk: &PublicKey,
        m: &Tn,
    ) -> Result<Self> {
        let tlev: Vec<TGLWE> = (1..l + 1)
            .map(|i| {
                TGLWE::encrypt(
                    &mut rng,
                    &param,
                    pk,
                    &(m * &(u64::MAX / beta.pow(i as u32) as u64)),
                )
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
        m: &Tn,
    ) -> Result<Self> {
        let tlev: Vec<TGLWE> = (1..l as u64 + 1)
            .map(|i| {
                let aux = if i < 64 {
                    m * &(u64::MAX / (1u64 << i))
                } else {
                    // 1<<64 would overflow, and anyways we're dividing u64::MAX
                    // by it, which would be equal to 1
                    m.clone() // TODO rm clone
                };
                TGLWE::encrypt_s(&mut rng, &param, sk, &aux)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }

    pub fn decrypt(&self, sk: &SecretKey, beta: u32) -> Tn {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, u64::MAX)
    }
}

impl TGLev {
    pub fn iter(&self) -> std::slice::Iter<TGLWE> {
        self.0.iter()
    }
}

// dot product between a TGLev and Vec<Tn<N>>, usually Vec<Tn<N>> comes from a
// decomposition of Tn<N>
// TGLev * Vec<Tn<N>> --> TGLWE
impl Mul<Vec<Tn>> for TGLev {
    type Output = TGLWE;
    fn mul(self, v: Vec<Tn>) -> Self::Output {
        assert_eq!(self.0.len(), v.len());

        // l TGLWES
        let tlwes: Vec<TGLWE> = self.0;
        let r: TGLWE = zip_eq(v, tlwes).map(|(a_d_i, glwe_i)| glwe_i * a_d_i).sum();
        r
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;

    use super::*;
    #[test]
    fn test_external_product() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 4,
            t: 16, // plaintext modulus
        };

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..50 {
            let (sk, _) = TGLWE::new_key(&mut rng, &param)?;

            let m1: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Tn = TGLev::encode(&param, &m1);

            let m2: Rq = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p2: Tn = TGLWE::encode(&param, &m2); // scaled by delta

            let tgsw = TGGSW::encrypt_s(&mut rng, &param, beta, l, &sk, &p1)?;
            let tlwe = TGLWE::encrypt_s(&mut rng, &param, &sk, &p2)?;

            let res: TGLWE = tgsw * tlwe;

            // let p_recovered = res.decrypt(&sk, beta);
            let p_recovered = res.decrypt(&sk);
            // downscaled by delta^-1
            let res_recovered = TGLWE::decode(&param, &p_recovered);

            // assert_eq!(m1 * m2, m_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), res_recovered);
        }

        Ok(())
    }
}
