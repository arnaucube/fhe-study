use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, RingParam, Rq, Tn, Zq, T64, TR};
use gfhe::{glwe, glwe::Param, GLWE};

use crate::tggsw::TGGSW;
use crate::tlev::TLev;
use crate::{tglwe, tglwe::TGLWE};

pub struct SecretKey(pub glwe::SecretKey<T64>);

impl SecretKey {
    /// from TFHE [2018-421] paper: A TLWE key k \in B^n, can be interpreted as a
    /// TRLWE key K \in B_N[X]^k having the same sequence of coefficients and
    /// vice-versa.
    pub fn to_tglwe(self, param: &Param) -> crate::tglwe::SecretKey {
        let s: TR<T64> = self.0 .0; // of length K*N
        assert_eq!(s.r.len(), param.k * param.ring.n); // sanity check

        // split into K vectors, and interpret each of them as a T_N[X]/(X^N+1)
        // polynomial
        let r: Vec<Tn> =
            s.r.chunks(param.ring.n)
                .map(|v| Tn::from_vec(&param.ring, v.to_vec()))
                .collect();
        crate::tglwe::SecretKey(glwe::SecretKey::<Tn>(TR { k: param.k, r }))
    }
}

pub type PublicKey = glwe::PublicKey<T64>;

#[derive(Clone, Debug)]
pub struct KSK(Vec<TLev>);

#[derive(Clone, Debug)]
pub struct TLWE(pub GLWE<T64>);

impl TLWE {
    pub fn zero(k: usize, ring_param: &RingParam) -> Self {
        Self(GLWE::<T64>::zero(k, ring_param))
    }

    pub fn new_key(rng: impl Rng, param: &Param) -> Result<(SecretKey, PublicKey)> {
        let (sk, pk): (glwe::SecretKey<T64>, glwe::PublicKey<T64>) = GLWE::new_key(rng, param)?;
        Ok((SecretKey(sk), pk))
    }

    pub fn encode(param: &Param, m: &Rq) -> T64 {
        assert_eq!(param.ring.n, 1);
        debug_assert_eq!(param.t, m.param.q); // plaintext modulus

        let delta = u64::MAX / param.t; // floored
        let coeffs = m.coeffs();
        T64(coeffs[0].v * delta)
    }
    pub fn decode(param: &Param, p: &T64) -> Rq {
        let p = p.mul_div_round(param.t, u64::MAX);
        Rq::from_vec_u64(&param.pt(), p.coeffs().iter().map(|c| c.0).collect())
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(rng: impl Rng, param: &Param, sk: &SecretKey, p: &T64) -> Result<Self> {
        let glwe = GLWE::encrypt_s(rng, param, &sk.0, p)?;
        Ok(Self(glwe))
    }
    pub fn encrypt(rng: impl Rng, param: &Param, pk: &PublicKey, p: &T64) -> Result<Self> {
        let glwe = GLWE::encrypt(rng, param, pk, p)?;
        Ok(Self(glwe))
    }
    pub fn decrypt(&self, sk: &SecretKey) -> T64 {
        self.0.decrypt(&sk.0)
    }

    pub fn new_ksk(
        mut rng: impl Rng,
        param: &Param,
        beta: u32,
        l: u32,
        sk: &SecretKey,
        new_sk: &SecretKey,
    ) -> Result<KSK> {
        let r: Vec<TLev> = (0..param.k)
            .into_iter()
            .map(|i|
                // treat sk_i as the msg being encrypted
                TLev::encrypt_s(&mut rng, param, beta, l, &new_sk, &sk.0.0 .r[i]))
            .collect::<Result<Vec<_>>>()?;

        Ok(KSK(r))
    }
    pub fn key_switch(&self, param: &Param, beta: u32, l: u32, ksk: &KSK) -> Self {
        let (a, b): (TR<T64>, T64) = (self.0 .0.clone(), self.0 .1);

        let lhs: TLWE = TLWE(GLWE(TR::zero(param.k * param.ring.n, &param.ring), b));

        // K iterations, ksk.0 contains K times GLev
        let rhs: TLWE = zip_eq(a.r, ksk.0.clone())
            .map(|(a_i, ksk_i)| ksk_i * a_i.decompose(beta, l)) // dot_product
            .sum();

        lhs - rhs
    }
    // modulus switch from Q (2^64) to Q2 (in blind_rotation Q2=K*N)
    pub fn mod_switch(&self, q2: u64) -> Self {
        let a: TR<T64> = self.0 .0.mod_switch(q2);
        let b: T64 = self.0 .1.mod_switch(q2);
        Self(GLWE(a, b))
    }
}

pub fn blind_rotation(
    param: &Param,
    c: TLWE, // kn
    btk: BootstrappingKey,
    table: TGLWE, // n,k
) -> TGLWE {
    debug_assert_eq!(c.0 .0.k, param.k);

    // TODO replace `param.k*param.ring.n` by `param.kn()`
    let c_kn: TLWE = c.mod_switch((param.k * param.ring.n) as u64);
    let (a, b): (TR<T64>, T64) = (c_kn.0 .0, c_kn.0 .1);
    // two main parts: rotate by a known power of X, rotate by a secret
    // power of X (using the C gate)

    // table * X^-b, ie. left rotate
    let v_xb: TGLWE = table.left_rotate(b.0 as usize);

    // rotate by a secret power of X using the cmux gate
    let mut c_j: TGLWE = v_xb.clone();
    let _ = (1..param.k).map(|j| {
        c_j = TGGSW::cmux(
            btk.0[j].clone(),
            c_j.clone(),
            c_j.clone().left_rotate(a.r[j].0 as usize),
        );
    });
    c_j
}

pub fn bootstrapping(
    param: &Param,
    btk: BootstrappingKey,
    table: TGLWE,
    c: TLWE, // kn
) -> TLWE {
    // kn
    let rotated: TGLWE = blind_rotation(param, c, btk.clone(), table);
    let c_h: TLWE = rotated.sample_extraction(&param, 0);
    let r = c_h.key_switch(param, 2, 64, &btk.1);
    r
}

#[derive(Clone, Debug)]
pub struct BootstrappingKey(
    pub Vec<TGGSW>,
    pub KSK, // kn
);
impl BootstrappingKey {
    pub fn from_sk(mut rng: impl Rng, param: &Param, sk: &tglwe::SecretKey) -> Result<Self> {
        let (beta, l) = (2u32, 64u32); // TMP

        let s: TR<Tn> = sk.0 .0.clone();
        let (sk2, _) = TLWE::new_key(&mut rng, &param.lwe())?; // TLWE<KN> compatible with TGLWE<N,K>

        // each btk_j = TGGSW_sk(s_i)
        let btk: Vec<TGGSW> = s
            .iter()
            .map(|s_i| TGGSW::encrypt_s(&mut rng, param, beta, l, sk, s_i))
            .collect::<Result<Vec<_>>>()?;

        let ksk = TLWE::new_ksk(
            &mut rng,
            &param.lwe(),
            beta,
            l,
            &sk.to_tlwe(&param.lwe()), // converted to length k*n
            &sk2,                      // created with length k*n
        )?;
        debug_assert_eq!(ksk.0.len(), param.lwe().k);
        debug_assert_eq!(ksk.0.len(), param.k * param.ring.n);

        Ok(Self(btk, ksk))
    }
}

pub fn compute_lookup_table(param: &Param) -> TGLWE {
    // from 2021-1402:
    // v(x) = \sum_j^{N-1} [(p_j / 2N mod p)/p] X^j

    // matrix of coefficients with size K*N = delta x T
    let delta: usize = param.ring.n / param.t as usize;
    let values: Vec<Zq> = (0..param.t).map(|v| Zq::from_u64(param.t, v)).collect();
    let coeffs: Vec<Zq> = (0..param.t as usize)
        .flat_map(|i| vec![values[i]; delta])
        .collect();
    let table = Rq::from_vec(&param.pt(), coeffs);

    // encode the table as plaintext
    let v: Tn = TGLWE::encode(param, &table);

    // encode the table as TGLWE ciphertext
    let v: TGLWE = TGLWE::from_plaintext(param.k, &param.ring, v);
    v
}

impl Add<TLWE> for TLWE {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        debug_assert_eq!(self.0 .0.k, other.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), other.0 .1.param());
        Self(self.0 + other.0)
    }
}
impl AddAssign for TLWE {
    fn add_assign(&mut self, rhs: Self) {
        debug_assert_eq!(self.0 .0.k, rhs.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), rhs.0 .1.param());
        self.0 += rhs.0
    }
}
impl Sum<TLWE> for TLWE {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let first = iter.next().unwrap();
        iter.fold(first, |acc, e| acc + e)
    }
}

impl Sub<TLWE> for TLWE {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        debug_assert_eq!(self.0 .0.k, other.0 .0.k);
        debug_assert_eq!(self.0 .1.param(), other.0 .1.param());
        Self(self.0 - other.0)
    }
}

// plaintext addition
impl Add<T64> for TLWE {
    type Output = Self;
    fn add(self, plaintext: T64) -> Self {
        let a: TR<T64> = self.0 .0;
        let b: T64 = self.0 .1 + plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext substraction
impl Sub<T64> for TLWE {
    type Output = Self;
    fn sub(self, plaintext: T64) -> Self {
        let a: TR<T64> = self.0 .0;
        let b: T64 = self.0 .1 - plaintext;
        Self(GLWE(a, b))
    }
}
// plaintext multiplication
impl Mul<T64> for TLWE {
    type Output = Self;
    fn mul(self, plaintext: T64) -> Self {
        let a: TR<T64> = TR {
            k: self.0 .0.k,
            r: self.0 .0.r.iter().map(|r_i| *r_i * plaintext).collect(),
        };
        let b: T64 = self.0 .1 * plaintext;
        Self(GLWE(a, b))
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use rand::distributions::Uniform;
    use std::time::Instant;

    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p: T64 = TLWE::encode(&param, &m);

            let c = TLWE::encrypt(&mut rng, &param, &pk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = TLWE::decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);

            // same but using encrypt_s (with sk instead of pk))
            let c = TLWE::encrypt_s(&mut rng, &param, &sk, &p)?;
            let p_recovered = c.decrypt(&sk);
            let m_recovered = TLWE::decode(&param, &p_recovered);

            assert_eq!(m, m_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLWE::encode(&param, &m1); // plaintext
            let p2: T64 = TLWE::encode(&param, &m2); // plaintext

            let c1 = TLWE::encrypt(&mut rng, &param, &pk, &p1)?;
            let c2 = TLWE::encrypt(&mut rng, &param, &pk, &p2)?;

            let c3 = c1 + c2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = TLWE::decode(&param, &p3_recovered);

            assert_eq!((m1 + m2).remodule(param.t), m3_recovered.remodule(param.t));
        }

        Ok(())
    }

    #[test]
    fn test_add_plaintext() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLWE::encode(&param, &m1); // plaintext
            let p2: T64 = TLWE::encode(&param, &m2); // plaintext

            let c1 = TLWE::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 + p2;

            let p3_recovered = c3.decrypt(&sk);
            let m3_recovered = TLWE::decode(&param, &p3_recovered);

            assert_eq!(m1 + m2, m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let mut rng = rand::thread_rng();
        let msg_dist = Uniform::new(0_u64, param.t);

        for _ in 0..200 {
            let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

            let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
            let p1: T64 = TLWE::encode(&param, &m1);
            // don't scale up p2, set it directly from m2
            let p2: T64 = T64(m2.coeffs()[0].v);

            let c1 = TLWE::encrypt(&mut rng, &param, &pk, &p1)?;

            let c3 = c1 * p2;

            let p3_recovered: T64 = c3.decrypt(&sk);
            let m3_recovered = TLWE::decode(&param, &p3_recovered);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq(param.t), m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_key_switch() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam { q: u64::MAX, n: 1 },
            k: 16,
            t: 128, // plaintext modulus
        };

        let beta: u32 = 2;
        let l: u32 = 64;

        let mut rng = rand::thread_rng();

        let (sk, pk) = TLWE::new_key(&mut rng, &param)?;
        let (sk2, _) = TLWE::new_key(&mut rng, &param)?;
        // ksk to switch from sk to sk2
        let ksk = TLWE::new_ksk(&mut rng, &param, beta, l, &sk, &sk2)?;

        let msg_dist = Uniform::new(0_u64, param.t);
        let m = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
        let p = TLWE::encode(&param, &m); // plaintext

        let c = TLWE::encrypt_s(&mut rng, &param, &sk, &p)?;

        let c2 = c.key_switch(&param, beta, l, &ksk);

        // decrypt with the 2nd secret key
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = TLWE::decode(&param, &p_recovered);
        assert_eq!(m.remodule(param.t), m_recovered.remodule(param.t));

        // do the same but now encrypting with pk
        let c = TLWE::encrypt(&mut rng, &param, &pk, &p)?;
        let c2 = c.key_switch(&param, beta, l, &ksk);
        let p_recovered = c2.decrypt(&sk2);
        let m_recovered = TLWE::decode(&param, &p_recovered);
        assert_eq!(m, m_recovered);

        Ok(())
    }

    #[test]
    fn test_bootstrapping() -> Result<()> {
        let param = Param {
            err_sigma: crate::ERR_SIGMA,
            ring: RingParam {
                q: u64::MAX,
                n: 1024,
            },
            k: 1,
            t: 128, // plaintext modulus
        };
        let mut rng = rand::thread_rng();

        let start = Instant::now();
        let table: TGLWE = compute_lookup_table(&param);
        println!("table took: {:?}", start.elapsed());

        let (sk, _) = TGLWE::new_key(&mut rng, &param)?;
        let sk_tlwe: SecretKey = sk.to_tlwe(&param);

        let start = Instant::now();
        let btk = BootstrappingKey::from_sk(&mut rng, &param, &sk)?;
        println!("btk took: {:?}", start.elapsed());

        let msg_dist = Uniform::new(0_u64, param.t);
        let m = Rq::rand_u64(&mut rng, msg_dist, &param.lwe().pt())?; // q=t, n=1
        let p = TLWE::encode(&param.lwe(), &m); // plaintext

        let c = TLWE::encrypt_s(&mut rng, &param.lwe(), &sk_tlwe, &p)?;

        let start = Instant::now();
        let bootstrapped: TLWE = bootstrapping(&param, btk, table, c);
        println!("bootstrapping took: {:?}", start.elapsed());

        let p_recovered: T64 = bootstrapped.decrypt(&sk_tlwe);
        let m_recovered = TLWE::decode(&param.lwe(), &p_recovered);
        assert_eq!(m_recovered, m);

        Ok(())
    }
}
