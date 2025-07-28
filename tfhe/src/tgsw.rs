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

    pub fn cmux(bit: Self, ct1: TLWE<K>, ct2: TLWE<K>) -> TLWE<K> {
        ct1.clone() + (bit * (ct2 - ct1))
    }
}

/// External product TGSW x TLWE
impl<const K: usize> Mul<TLWE<K>> for TGSW<K> {
    type Output = TLWE<K>;

    fn mul(self, tlwe: TLWE<K>) -> TLWE<K> {
        let beta: u32 = 2;
        let l: u32 = 64; // TODO wip

        // since N=1, each tlwe element is a vector of length=1, decomposed into
        // l elements, and we have K of them
        let tlwe_ab: Vec<T64> = [tlwe.0 .0 .0.clone(), vec![tlwe.0 .1]].concat();

        let tgsw_ab: Vec<TLev<K>> = [self.0.clone(), vec![self.1]].concat();
        assert_eq!(tgsw_ab.len(), tlwe_ab.len());

        let r: TLWE<K> = zip_eq(tgsw_ab, tlwe_ab)
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

    #[test]
    fn test_external_product() -> Result<()> {
        const T: u64 = 2; // plaintext modulus
        const K: usize = 32;

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..50 {
            let (sk, _) = TLWE::<K>::new_key(&mut rng)?;

            let m1: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p1: T64 = TLev::<K>::encode::<T>(&m1);

            let m2: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p2: T64 = TLWE::<K>::encode::<T>(&m2); // scaled by delta

            let tgsw = TGSW::<K>::encrypt_s(&mut rng, beta, l, &sk, &p1)?;
            let tlwe = TLWE::<K>::encrypt_s(&mut rng, &sk, &p2)?;

            let res: TLWE<K> = tgsw * tlwe;

            // let p_recovered = res.decrypt(&sk, beta);
            let p_recovered = res.decrypt(&sk);
            // downscaled by delta^-1
            let res_recovered = TLWE::<K>::decode::<T>(&p_recovered);

            // assert_eq!(m1 * m2, m_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), res_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_cmux() -> Result<()> {
        const T: u64 = 2; // plaintext modulus
        const K: usize = 32;

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, T);

        for _ in 0..50 {
            let (sk, _) = TLWE::<K>::new_key(&mut rng)?;

            let m1: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p1: T64 = TLWE::<K>::encode::<T>(&m1); // scaled by delta

            let m2: Rq<T, 1> = Rq::rand_u64(&mut rng, msg_dist)?;
            let p2: T64 = TLWE::<K>::encode::<T>(&m2); // scaled by delta

            for bit_raw in 0..2 {
                let bit = TGSW::<K>::encrypt_s(&mut rng, beta, l, &sk, &T64(bit_raw))?;

                let c1 = TLWE::<K>::encrypt_s(&mut rng, &sk, &p1)?;
                let c2 = TLWE::<K>::encrypt_s(&mut rng, &sk, &p2)?;

                let res: TLWE<K> = TGSW::cmux(bit, c1, c2);

                let p_recovered = res.decrypt(&sk);
                // downscaled by delta^-1
                let res_recovered = TLWE::<K>::decode::<T>(&p_recovered);

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
