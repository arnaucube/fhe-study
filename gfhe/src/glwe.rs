use anyhow::Result;
use itertools::zip_eq;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub};

use arith::{Ring, Rq, Zq, TR};

use crate::glev::GLev;

const ERR_SIGMA: f64 = 3.2;

#[derive(Clone, Debug)]
pub struct GLWE<const Q: u64, const N: usize, const K: usize>(TR<Rq<Q, N>, K>, Rq<Q, N>);

#[derive(Clone, Debug)]
pub struct SecretKey<const Q: u64, const N: usize, const K: usize>(TR<Rq<Q, N>, K>);
#[derive(Clone, Debug)]
pub struct PublicKey<const Q: u64, const N: usize, const K: usize>(Rq<Q, N>, TR<Rq<Q, N>, K>);

// K GLevs, each KSK_i=l GLWEs
#[derive(Clone, Debug)]
pub struct KSK<const Q: u64, const N: usize, const K: usize>(Vec<GLev<Q, N, K>>);

impl<const Q: u64, const N: usize, const K: usize> GLWE<Q, N, K> {
    pub fn zero() -> Self {
        Self(TR::zero(), Rq::zero())
    }

    pub fn new_key(mut rng: impl Rng) -> Result<(SecretKey<Q, N, K>, PublicKey<Q, N, K>)> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let s: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Xi_key);
        let a: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Uniform::new(0_f64, Q as f64));
        let e = Rq::<Q, N>::rand(&mut rng, Xi_err);

        let pk: PublicKey<Q, N, K> = PublicKey((&a * &s) + e, a);
        Ok((SecretKey(s), pk))
    }

    pub fn new_ksk(
        mut rng: impl Rng,
        beta: u32,
        l: u32,
        sk: &SecretKey<Q, N, K>,
        new_sk: &SecretKey<Q, N, K>,
    ) -> Result<KSK<Q, N, K>> {
        let r: Vec<GLev<Q, N, K>> = (0..K)
            .into_iter()
            .map(|i|
                // treat sk_i as the msg being encrypted
                GLev::<Q, N, K>::encrypt_s(&mut rng, beta, l, &new_sk, &sk.0 .0[i]))
            .collect::<Result<Vec<_>>>()?;

        Ok(KSK(r))
    }
    pub fn key_switch(&self, beta: u32, l: u32, ksk: &KSK<Q, N, K>) -> Self {
        let (a, b): (TR<Rq<Q, N>, K>, Rq<Q, N>) = (self.0.clone(), self.1);

        let lhs: GLWE<Q, N, K> = GLWE(TR::zero(), b);

        // K iterations, ksk.0 contains K times GLev
        let rhs: GLWE<Q, N, K> = zip_eq(a.0, ksk.0.clone())
            .map(|(a_i, ksk_i)| Self::dot_prod(a_i.decompose(beta, l), ksk_i))
            .sum();

        lhs - rhs
    }
    // note: a_decomp is of length N
    fn dot_prod(a_decomp: Vec<Rq<Q, N>>, ksk_i: GLev<Q, N, K>) -> GLWE<Q, N, K> {
        // l times GLWES
        let glwes: Vec<GLWE<Q, N, K>> = ksk_i.0;

        // l iterations
        let r: GLWE<Q, N, K> = zip_eq(a_decomp, glwes)
            .map(|(a_d_i, glwe_i)| glwe_i * a_d_i)
            .sum();
        r
    }

    // encrypts with the given SecretKey (instead of PublicKey)
    pub fn encrypt_s(
        mut rng: impl Rng,
        sk: &SecretKey<Q, N, K>,
        m: &Rq<Q, N>,
        // TODO delta not as input
        delta: u64,
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Xi_key);
        let e = Rq::<Q, N>::rand(&mut rng, Xi_err);

        let b: Rq<Q, N> = (&a * &sk.0) + *m * delta + e;
        Ok(Self(a, b))
    }
    pub fn encrypt(
        mut rng: impl Rng,
        pk: &PublicKey<Q, N, K>,
        m: &Rq<Q, N>,
        delta: u64,
    ) -> Result<Self> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u: Rq<Q, N> = Rq::rand(&mut rng, Xi_key);

        let e0 = Rq::<Q, N>::rand(&mut rng, Xi_err);
        let e1 = TR::<Rq<Q, N>, K>::rand(&mut rng, Xi_err);

        let b: Rq<Q, N> = pk.0 * u + *m * delta + e0;
        let d: TR<Rq<Q, N>, K> = &pk.1 * &u + e1;

        Ok(Self(d, b))
    }
    pub fn decrypt<const T: u64>(&self, sk: &SecretKey<Q, N, K>, delta: u64) -> Rq<Q, N> {
        let (d, b): (TR<Rq<Q, N>, K>, Rq<Q, N>) = (self.0.clone(), self.1);
        let r: Rq<Q, N> = b - &d * &sk.0;
        let r = r.mul_div_round(T, Q);
        r
    }

    pub fn mod_switch<const P: u64>(&self) -> GLWE<P, N, K> {
        let a: TR<Rq<P, N>, K> = TR(self.0 .0.iter().map(|r| r.mod_switch::<P>()).collect());
        let b: Rq<P, N> = self.1.mod_switch::<P>();
        GLWE(a, b)
    }
}

impl<const Q: u64, const N: usize, const K: usize> Add<GLWE<Q, N, K>> for GLWE<Q, N, K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let a: TR<Rq<Q, N>, K> = self.0 + other.0;
        let b: Rq<Q, N> = self.1 + other.1;
        Self(a, b)
    }
}

impl<const Q: u64, const N: usize, const K: usize> Add<Rq<Q, N>> for GLWE<Q, N, K> {
    type Output = Self;
    fn add(self, plaintext: Rq<Q, N>) -> Self {
        let a: TR<Rq<Q, N>, K> = self.0;
        let b: Rq<Q, N> = self.1 + plaintext;
        Self(a, b)
    }
}
impl<const Q: u64, const N: usize, const K: usize> AddAssign for GLWE<Q, N, K> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..K {
            self.0 .0[i] = self.0 .0[i] + rhs.0 .0[i];
        }
        self.1 = self.1 + rhs.1;
    }
}
impl<const Q: u64, const N: usize, const K: usize> Sum<GLWE<Q, N, K>> for GLWE<Q, N, K> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = GLWE::<Q, N, K>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<const Q: u64, const N: usize, const K: usize> Sub<GLWE<Q, N, K>> for GLWE<Q, N, K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let a: TR<Rq<Q, N>, K> = self.0 - other.0;
        let b: Rq<Q, N> = self.1 - other.1;
        Self(a, b)
    }
}

impl<const Q: u64, const N: usize, const K: usize> Mul<Rq<Q, N>> for GLWE<Q, N, K> {
    type Output = Self;
    fn mul(self, plaintext: Rq<Q, N>) -> Self {
        // first compute the NTT for plaintext, to avoid computing it at each
        // iteration, speeding up the multiplications
        let mut plaintext = plaintext.clone();
        plaintext.compute_evals();

        let a: TR<Rq<Q, N>, K> = TR(self.0 .0.iter().map(|r_i| *r_i * plaintext).collect());
        let b: Rq<Q, N> = self.1 * plaintext;
        Self(a, b)
    }
}

impl<const Q: u64, const N: usize, const K: usize> Mul<Zq<Q>> for GLWE<Q, N, K> {
    type Output = Self;
    fn mul(self, e: Zq<Q>) -> Self {
        let a: TR<Rq<Q, N>, K> = TR(self.0 .0.iter().map(|r_i| *r_i * e).collect());
        let b: Rq<Q, N> = self.1 * e;
        Self(a, b)
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
        const T: u64 = 32; // plaintext modulus
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m: Rq<Q, N> = m.remodule::<Q>();

            let c = S::encrypt(&mut rng, &pk, &m, delta)?;
            let m_recovered = c.decrypt::<T>(&sk, delta);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());

            // same but using encrypt_s (with sk instead of pk))
            let c = S::encrypt_s(&mut rng, &sk, &m, delta)?;
            let m_recovered = c.decrypt::<T>(&sk, delta);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_addition() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 20;
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m1: Rq<Q, N> = m1.remodule::<Q>();
            let m2: Rq<Q, N> = m2.remodule::<Q>();

            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;
            let c2 = S::encrypt(&mut rng, &pk, &m2, delta)?;

            let c3 = c1 + c2;

            let m3_recovered = c3.decrypt::<T>(&sk, delta);

            assert_eq!((m1 + m2).remodule::<T>(), m3_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_add_plaintext() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 32;
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m1: Rq<Q, N> = m1.remodule::<Q>();
            let m2: Rq<Q, N> = m2.remodule::<Q>();
            let m2_scaled: Rq<Q, N> = m2 * delta;

            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;

            let c3 = c1 + m2_scaled;

            let m3_recovered = c3.decrypt::<T>(&sk, delta);

            assert_eq!((m1 + m2).remodule::<T>(), m3_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_mul_plaintext() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 16;
        const T: u64 = 4;
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m1 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m2 = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m1: Rq<Q, N> = m1.remodule::<Q>();
            let m2: Rq<Q, N> = m2.remodule::<Q>();
            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;

            let c3 = c1 * m2;

            let m3_recovered: Rq<Q, N> = c3.decrypt::<T>(&sk, delta);
            let m3_recovered: Rq<T, N> = m3_recovered.remodule::<T>();
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), m3_recovered);
        }

        Ok(())
    }

    #[test]
    fn test_mod_switch() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const P: u64 = 2u64.pow(8) + 1;
        // note: wip, Q and P chosen so that P/Q is an integer
        const N: usize = 8;
        const T: u64 = 8; // plaintext modulus, must be a prime or power of a prime
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        dbg!(P as f64 / Q as f64);
        dbg!(delta);
        dbg!(delta as f64 * P as f64 / Q as f64);
        dbg!(delta as f64 * (P as f64 / Q as f64));

        for _ in 0..200 {
            let (sk, pk) = S::new_key(&mut rng)?;

            let msg_dist = Uniform::new(0_u64, T);
            let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
            let m: Rq<Q, N> = m.remodule::<Q>();

            let c = S::encrypt(&mut rng, &pk, &m, delta)?;
            // let c = S::encrypt_s(&mut rng, &sk, &m, delta)?;

            let c2 = c.mod_switch::<P>();
            let sk2: SecretKey<P, N, K> =
                SecretKey(TR(sk.0 .0.iter().map(|s_i| s_i.remodule::<P>()).collect()));
            let delta2: u64 = ((P as f64 * delta as f64) / Q as f64).round() as u64;

            let m_recovered = c2.decrypt::<T>(&sk2, delta2);

            assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());
        }

        Ok(())
    }

    #[test]
    fn test_key_switch() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 128;
        const T: u64 = 2; // plaintext modulus
        const K: usize = 16;
        type S = GLWE<Q, N, K>;

        let beta: u32 = 2;
        let l: u32 = 16;

        let delta: u64 = Q / T; // floored
        let mut rng = rand::thread_rng();

        let (sk, pk) = S::new_key(&mut rng)?;
        let (sk2, _) = S::new_key(&mut rng)?;
        // ksk to switch from sk to sk2
        let ksk = S::new_ksk(&mut rng, beta, l, &sk, &sk2)?;

        let msg_dist = Uniform::new(0_u64, T);
        let m = Rq::<T, N>::rand_u64(&mut rng, msg_dist)?;
        let m: Rq<Q, N> = m.remodule::<Q>();

        let c = S::encrypt_s(&mut rng, &sk, &m, delta)?;

        let c2 = c.key_switch(beta, l, &ksk);

        // decrypt with the 2nd secret key
        let m_recovered = c2.decrypt::<T>(&sk2, delta);
        assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());

        // do the same but now encrypting with pk
        // let c = S::encrypt(&mut rng, &pk, &m, delta)?;
        // let c2 = c.key_switch(beta, l, &ksk);
        // let m_recovered = c2.decrypt::<T>(&sk2, delta);
        // assert_eq!(m.remodule::<T>(), m_recovered.remodule::<T>());

        Ok(())
    }
}
