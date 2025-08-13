use anyhow::Result;
use itertools::zip_eq;
use rand::distributions::Standard;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::array;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, RingParam, Rq, Tn, T64, TR};
use gfhe::{glwe, glwe::Param, GLWE};

use crate::tlev::TLev;
use crate::{tlwe, tlwe::TLWE};

// pub type SecretKey<const N: usize, const K: usize> = glwe::SecretKey<Tn<N>, K>;
#[derive(Clone, Debug)]
pub struct SecretKey(pub glwe::SecretKey<Tn>);
// pub struct SecretKey<const K: usize>(pub tlwe::SecretKey<K>);

impl SecretKey {
    pub fn to_tlwe(&self, param: &Param) -> tlwe::SecretKey {
        let s: TR<Tn> = self.0 .0.clone();
        debug_assert_eq!(s.r.len(), param.k); // sanity check

        let kn = param.k * param.ring.n;
        let r: Vec<Vec<T64>> = s.r.iter().map(|s_i| s_i.coeffs()).collect();
        let r: Vec<T64> = r.into_iter().flatten().collect();
        debug_assert_eq!(r.len(), kn); // sanity check
        tlwe::SecretKey(glwe::SecretKey::<T64>(TR::<T64>::new(kn, r)))
    }
}

pub type PublicKey = glwe::PublicKey<Tn>;

#[derive(Clone, Debug)]
pub struct TGLWE(pub GLWE<Tn>);

impl TGLWE {
    pub fn zero(k: usize, param: &RingParam) -> Self {
        Self(GLWE::<Tn>::zero(k, param))
    }
    pub fn from_plaintext(k: usize, param: &RingParam, p: Tn) -> Self {
        Self(GLWE::<Tn>::from_plaintext(k, param, p))
    }

    pub fn new_key(mut rng: impl Rng, param: &Param) -> Result<(SecretKey, PublicKey)> {
        // assert_eq!(KN, K * N); // this is wip, while not being able to compute K*N
        let (sk_tlwe, _) = TLWE::new_key(&mut rng, &param.lwe())?; //param.lwe() so that it uses K*N
        debug_assert_eq!(sk_tlwe.0 .0.r.len(), param.lwe().k); // =KN (sanity check)

        // let sk = crate::tlwe::sk_to_tglwe::<N, K, KN>(sk_tlwe);
        let sk = sk_tlwe.to_tglwe(param);
        let pk: PublicKey = GLWE::pk_from_sk(rng, param, sk.0.clone())?;
        Ok((sk, pk))
    }

    pub fn encode(param: &Param, m: &Rq) -> Tn {
        debug_assert_eq!(param.t, m.param.q); // plaintext modulus
        let p = param.t;
        let delta = u64::MAX / p; // floored
        let coeffs = m.coeffs();
        // Tn(array::from_fn(|i| T64(coeffs[i].0 * delta)))
        Tn {
            param: param.ring,
            coeffs: coeffs.iter().map(|c_i| T64(c_i.v * delta)).collect(),
        }
    }
    pub fn decode(param: &Param, pt: &Tn) -> Rq {
        let p = param.t;
        let pt = pt.mul_div_round(p, u64::MAX);
        Rq::from_vec_u64(&param.pt(), pt.coeffs().iter().map(|c| c.0).collect())
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(rng: impl Rng, param: &Param, sk: &SecretKey, p: &Tn) -> Result<Self> {
        let glwe = GLWE::encrypt_s(rng, param, &sk.0, p)?;
        Ok(Self(glwe))
    }
    pub fn encrypt(rng: impl Rng, param: &Param, pk: &PublicKey, p: &Tn) -> Result<Self> {
        let glwe = GLWE::encrypt(rng, param, &pk, p)?;
        Ok(Self(glwe))
    }
    pub fn decrypt(&self, sk: &SecretKey) -> Tn {
        self.0.decrypt(&sk.0)
    }

    /// Sample extraction / Coefficient extraction
    pub fn sample_extraction(&self, param: &Param, h: usize) -> TLWE {
        let n = param.ring.n;
        assert!(h < n);

        let a: TR<Tn> = self.0 .0.clone();
        // set a_{n*i+j} = a_{i, h-j}     if j \in {0, h}
        //                 -a_{i, n+h-j}  if j \in {h+1, n-1}
        let new_a: Vec<T64> = a
            .iter()
            .flat_map(|a_i| {
                let a_i = a_i.coeffs();
                (0..n)
                    .map(|j| if j <= h { a_i[h - j] } else { -a_i[n + h - j] })
                    .collect::<Vec<T64>>()
            })
            .collect::<Vec<T64>>();
        debug_assert_eq!(new_a.len(), param.k * param.ring.n); // sanity check

        TLWE(GLWE(
            TR {
                // TODO use constructor `new`, which will check len with k
                k: param.k * param.ring.n,
                r: new_a,
            },
            self.0 .1.coeffs()[h],
        ))
    }
    pub fn left_rotate(&self, h: usize) -> Self {
        let (a, b): (TR<Tn>, Tn) = (self.0 .0.clone(), self.0 .1.clone());
        Self(GLWE(a.left_rotate(h), b.left_rotate(h)))
    }
}

impl Add<TGLWE> for TGLWE {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        debug_assert_eq!(self.0 .0.k, other.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), other.0 .1.param());

        Self(self.0 + other.0)
    }
}
impl AddAssign for TGLWE {
    fn add_assign(&mut self, other: Self) {
        debug_assert_eq!(self.0 .0.k, other.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), other.0 .1.param());

        self.0 += other.0
    }
}
impl Sum<TGLWE> for TGLWE {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        // let mut acc = TGLWE::<N, K>::zero();
        // for e in iter {
        //     acc += e;
        // }
        // acc
        let first = iter.next().unwrap();
        iter.fold(first, |acc, e| acc + e)
    }
}

impl Sub<TGLWE> for TGLWE {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        debug_assert_eq!(self.0 .0.k, other.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), other.0 .1.param());

        Self(self.0 - other.0)
    }
}

// plaintext addition
impl Add<Tn> for TGLWE {
    type Output = Self;
    fn add(self, plaintext: Tn) -> Self {
        debug_assert_eq!(self.0 .1.param(), plaintext.param());

        let a: TR<Tn> = self.0 .0;
        let b: Tn = self.0 .1 + plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext substraction
impl Sub<Tn> for TGLWE {
    type Output = Self;
    fn sub(self, plaintext: Tn) -> Self {
        debug_assert_eq!(self.0 .1.param(), plaintext.param());

        let a: TR<Tn> = self.0 .0;
        let b: Tn = self.0 .1 - plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext multiplication
impl Mul<Tn> for TGLWE {
    type Output = Self;
    fn mul(self, plaintext: Tn) -> Self {
        debug_assert_eq!(self.0 .1.param(), plaintext.param());

        let a: TR<Tn> = TR {
            k: self.0 .0.k,
            r: self.0 .0.r.iter().map(|r_i| r_i * &plaintext).collect(),
        };
        let b: Tn = self.0 .1 * plaintext;
        Self(GLWE(a, b))
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
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TGLWE::new_key(&mut rng, &param)?;

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p: Tn = TGLWE::encode(&param, &m);

            let c = TGLWE::encrypt(&mut rng, &param, &pk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = TGLWE::decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = TGLWE::encrypt_s(&mut rng, &param, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = TGLWE::decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TGLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Tn = TGLWE::encode(&param, &m1); // plaintext
            let p2: Tn = TGLWE::encode(&param, &m2); // plaintext

            let c1 = TGLWE::encrypt(&mut rng, &param, &pk, &p1)?;
            let c2 = TGLWE::encrypt(&mut rng, &param, &pk, &p2)?;

            let c3 = c1 + c2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = TGLWE::decode(&param, &p3_recovered);

            assert_eq!((m1 + m2).remodule(param.t), m3_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_add_plaintext() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TGLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Tn = TGLWE::encode(&param, &m1); // plaintext
            let p2: Tn = TGLWE::encode(&param, &m2); // plaintext

            let c1 = TGLWE::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = TGLWE::decode(&param, &p3_recovered);

            assert_eq!(m1 + m2, m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TGLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: Tn = TGLWE::encode(&param, &m1);
            // don't scale up p2, set it directly from m2
            let p2: Tn = Tn {
                param: param.ring,
                coeffs: m2.coeffs().iter().map(|c_i| T64(c_i.v)).collect(),
            };

            let c1 = TGLWE::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: Tn = c3.decrypt(&sk);
            let m3_recovered = TGLWE::decode(&param, &p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_sample_extraction() -> Result<()> {
        let param = Param {
            ring: RingParam { q: u64::MAX, n: 64 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..20 {
            let (sk, pk) = TGLWE::new_key(&mut rng, &param)?;
            let sk_tlwe = sk.to_tlwe(&param);

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p: Tn = TGLWE::encode(&param, &m);

            let c = TGLWE::encrypt(&mut rng, &param, &pk, &p)?;

            for h in 0..param.ring.n {
                let c_h: TLWE = c.sample_extraction(&param, h);

                let p_recovered = c_h.decrypt(&sk_tlwe);
                let m_recovered = TLWE::decode(&param, &p_recovered);
                assert_eq!(m.coeffs()[h], m_recovered.coeffs()[0]);
            }
        }

        Ok(())
    }
}
