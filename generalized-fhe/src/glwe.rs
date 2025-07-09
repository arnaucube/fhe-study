use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};
use std::ops::{Add, Mul};

use arith::{Ring, Rq, TR};

const ERR_SIGMA: f64 = 3.2;

pub struct GLWE<const Q: u64, const N: usize, const K: usize>(TR<Rq<Q, N>, K>, Rq<Q, N>);

#[derive(Clone, Debug)]
pub struct SecretKey<const Q: u64, const N: usize, const K: usize>(TR<Rq<Q, N>, K>);
#[derive(Clone, Debug)]
pub struct PublicKey<const Q: u64, const N: usize, const K: usize>(Rq<Q, N>, TR<Rq<Q, N>, K>);

impl<const Q: u64, const N: usize, const K: usize> GLWE<Q, N, K> {
    pub fn new_key(mut rng: impl Rng) -> Result<(SecretKey<Q, N, K>, PublicKey<Q, N, K>)> {
        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let s: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Xi_key);
        let a: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Uniform::new(0_f64, Q as f64));
        let e = Rq::<Q, N>::rand(&mut rng, Xi_err);

        let pk: PublicKey<Q, N, K> = PublicKey((&a * &s) + e, a);
        Ok((SecretKey(s), pk))
    }

    // TODO delta not as input
    pub fn encrypt_s<const T: u64>(
        mut rng: impl Rng,
        sk: &SecretKey<Q, N, K>,
        m: &Rq<T, N>,
        delta: u64,
    ) -> Result<Self> {
        let m: Rq<Q, N> = m.remodule::<Q>();

        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let a: TR<Rq<Q, N>, K> = TR::rand(&mut rng, Xi_key);
        let e = Rq::<Q, N>::rand(&mut rng, Xi_err);

        let b: Rq<Q, N> = (&a * &sk.0) + m * delta + e;
        Ok(Self(a, b))
    }
    pub fn encrypt<const T: u64>(
        mut rng: impl Rng,
        pk: &PublicKey<Q, N, K>,
        m: &Rq<T, N>,
        delta: u64,
    ) -> Result<Self> {
        let m: Rq<Q, N> = m.remodule::<Q>();

        let Xi_key = Uniform::new(0_f64, 2_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let u: Rq<Q, N> = Rq::rand(&mut rng, Xi_key);

        let e0 = Rq::<Q, N>::rand(&mut rng, Xi_err);
        let e1 = TR::<Rq<Q, N>, K>::rand(&mut rng, Xi_err);

        let b: Rq<Q, N> = pk.0 * u + m * delta + e0;
        let d: TR<Rq<Q, N>, K> = &pk.1 * &u + e1;

        Ok(Self(d, b))
    }
    pub fn decrypt<const T: u64>(&self, sk: &SecretKey<Q, N, K>, delta: u64) -> Rq<T, N> {
        let (d, b): (TR<Rq<Q, N>, K>, Rq<Q, N>) = (self.0.clone(), self.1);
        let r: Rq<Q, N> = b - &d * &sk.0;
        let r_scaled: Vec<f64> = r
            .coeffs()
            .iter()
            .map(|e| (e.0 as f64 / delta as f64).round())
            .collect();
        let r = Rq::<Q, N>::from_vec_f64(r_scaled);
        r.remodule::<T>()
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

            let c = S::encrypt(&mut rng, &pk, &m, delta)?;
            let m_recovered = c.decrypt(&sk, delta);

            assert_eq!(m, m_recovered);
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

            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;
            let c2 = S::encrypt(&mut rng, &pk, &m2, delta)?;

            let c3 = c1 + c2;

            let m3_recovered = c3.decrypt(&sk, delta);

            assert_eq!(m1 + m2, m3_recovered);
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
            let m2_scaled: Rq<Q, N> = m2.remodule::<Q>() * delta;

            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;

            let c3 = c1 + m2_scaled;

            let m3_recovered = c3.decrypt(&sk, delta);

            assert_eq!(m1 + m2, m3_recovered);
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
            let m2: Rq<Q, N> = m2.remodule::<Q>();
            let c1 = S::encrypt(&mut rng, &pk, &m1, delta)?;

            let c3 = c1 * m2;

            let m3_recovered: Rq<T, N> = c3.decrypt(&sk, delta);
            assert_eq!((m1.to_r() * m2.to_r()).to_rq::<T>(), m3_recovered);
        }

        Ok(())
    }
}
