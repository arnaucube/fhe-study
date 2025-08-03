use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::array;
use std::ops::{Add, Mul};

use arith::{Ring, Rq, Tn, T64, TR};

use crate::tglwe::{PublicKey, SecretKey, TGLWE};
use gfhe::glwe::GLWE;

/// vector of length K+1 = ([K * TGLev], [1 * TGLev])
#[derive(Clone, Debug)]
pub struct TGGSW<const N: usize, const K: usize>(pub(crate) Vec<TGLev<N, K>>, TGLev<N, K>);

impl<const N: usize, const K: usize> TGGSW<N, K> {
    pub fn encrypt_s(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<N, K>,
        m: &Tn<N>,
    ) -> Result<Self> {
        let a: Vec<TGLev<N, K>> = (0..K)
            .map(|i| TGLev::encrypt_s(&mut rng, beta, l, sk, &(-sk.0 .0 .0[i] * *m)))
            .collect::<Result<Vec<_>>>()?;
        let b: TGLev<N, K> = TGLev::encrypt_s(&mut rng, beta, l, sk, m)?;
        Ok(Self(a, b))
    }

    pub fn decrypt(&self, sk: &SecretKey<N, K>, beta: u32) -> Tn<N> {
        self.1.decrypt(sk, beta)
    }

    pub fn cmux(bit: Self, ct1: TGLWE<N, K>, ct2: TGLWE<N, K>) -> TGLWE<N, K> {
        ct1.clone() + (bit * (ct2 - ct1))
    }
}

/// External product TGGSW x TGLWE
impl<const N: usize, const K: usize> Mul<TGLWE<N, K>> for TGGSW<N, K> {
    type Output = TGLWE<N, K>;

    fn mul(self, tglwe: TGLWE<N, K>) -> TGLWE<N, K> {
        let beta: u32 = 2;
        let l: u32 = 64; // TODO wip

        let tglwe_ab: Vec<Tn<N>> = [tglwe.0 .0 .0.clone(), vec![tglwe.0 .1]].concat();

        let tgsw_ab: Vec<TGLev<N, K>> = [self.0.clone(), vec![self.1]].concat();
        assert_eq!(tgsw_ab.len(), tglwe_ab.len());

        let r: TGLWE<N, K> = zip_eq(tgsw_ab, tglwe_ab)
            .map(|(tlev_i, tglwe_i)| tlev_i * tglwe_i.decompose(beta, l))
            .sum();
        r
    }
}

#[derive(Clone, Debug)]
pub struct TGLev<const N: usize, const K: usize>(pub(crate) Vec<TGLWE<N, K>>);

impl<const N: usize, const K: usize> TGLev<N, K> {
    pub fn encode<const T: u64>(m: &Rq<T, N>) -> Tn<N> {
        let coeffs = m.coeffs();
        Tn(array::from_fn(|i| T64(coeffs[i].0)))
    }
    pub fn decode<const T: u64>(p: &Tn<N>) -> Rq<T, N> {
        Rq::<T, N>::from_vec_u64(p.coeffs().iter().map(|c| c.0).collect())
    }
    pub fn encrypt(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        pk: &PublicKey<N, K>,
        m: &Tn<N>,
    ) -> Result<Self> {
        let tlev: Vec<TGLWE<N, K>> = (1..l + 1)
            .map(|i| {
                TGLWE::<N, K>::encrypt(&mut rng, pk, &(*m * (u64::MAX / beta.pow(i as u32) as u64)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }
    pub fn encrypt_s(
        mut rng: impl Rng,
        _beta: u32, // TODO rm, and make beta=2 always
        l: u32,
        sk: &SecretKey<N, K>,
        m: &Tn<N>,
    ) -> Result<Self> {
        let tlev: Vec<TGLWE<N, K>> = (1..l as u64 + 1)
            .map(|i| {
                let aux = if i < 64 {
                    *m * (u64::MAX / (1u64 << i))
                } else {
                    // 1<<64 would overflow, and anyways we're dividing u64::MAX
                    // by it, which would be equal to 1
                    *m
                };
                TGLWE::<N, K>::encrypt_s(&mut rng, sk, &aux)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(tlev))
    }

    pub fn decrypt(&self, sk: &SecretKey<N, K>, beta: u32) -> Tn<N> {
        let pt = self.0[0].decrypt(sk);
        pt.mul_div_round(beta as u64, u64::MAX)
    }
}

impl<const N: usize, const K: usize> TGLev<N, K> {
    pub fn iter(&self) -> std::slice::Iter<TGLWE<N, K>> {
        self.0.iter()
    }
}

// dot product between a TGLev and Vec<Tn<N>>, usually Vec<Tn<N>> comes from a
// decomposition of Tn<N>
// TGLev * Vec<Tn<N>> --> TGLWE
impl<const N: usize, const K: usize> Mul<Vec<Tn<N>>> for TGLev<N, K> {
    type Output = TGLWE<N, K>;
    fn mul(self, v: Vec<Tn<N>>) -> Self::Output {
        assert_eq!(self.0.len(), v.len());

        // l TGLWES
        let tlwes: Vec<TGLWE<N, K>> = self.0;
        let r: TGLWE<N, K> = zip_eq(v, tlwes).map(|(a_d_i, glwe_i)| glwe_i * a_d_i).sum();
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
        const T: u64 = 16; // plaintext modulus
        const K: usize = 4;
        const N: usize = 64;
        const KN: usize = K * N;

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..50 {
            let (sk, _) = TGLWE::<N, K>::new_key::<KN>(&mut rng)?;

            let m1: Rq<T, N> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p1: Tn<N> = TGLev::<N, K>::encode::<T>(&m1);

            let m2: Rq<T, N> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p2: Tn<N> = TGLWE::<N, K>::encode::<T>(&m2); // scaled by delta

            let tgsw = TGGSW::<N, K>::encrypt_s(&mut rng, beta, l, &sk, &p1)?;
            let tlwe = TGLWE::<N, K>::encrypt_s(&mut rng, &sk, &p2)?;

            let res: TGLWE<N, K> = tgsw * tlwe;

            // let p_recovered = res.decrypt(&sk, beta);
            let p_recovered = res.decrypt(&sk);
            // downscaled by delta^-1
            let res_recovered = TGLWE::<N, K>::decode::<T>(&p_recovered);

            // assert_eq!(m1 * m2, m_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), res_recovered);
        }

        Ok(())
    }
}
