//! Implementation of CKKS https://eprint.iacr.org/2016/421.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

use arith::{RingParam, Rq, C, R};

use anyhow::Result;
use rand::Rng;
use rand_distr::{Normal, Uniform};

pub mod encoder;
pub use encoder::Encoder;

// error deviation for the Gaussian(Normal) distribution
// sigma=3.2 from: https://eprint.iacr.org/2016/421.pdf page 17
const ERR_SIGMA: f64 = 3.2;

#[derive(Clone, Copy, Debug)]
pub struct Params {
    ring: RingParam,
    t: u64,
}

#[derive(Debug)]
pub struct PublicKey(Rq, Rq);

pub struct SecretKey(Rq);

pub struct CKKS {
    params: Params,
    encoder: Encoder,
}

impl CKKS {
    pub fn new(params: &Params, delta: C<f64>) -> Self {
        let encoder = Encoder::new(params.ring.n, delta);
        Self {
            params: params.clone(),
            encoder,
        }
    }
    /// generate a new key pair (privK, pubK)
    pub fn new_key(&self, mut rng: impl Rng) -> Result<(SecretKey, PublicKey)> {
        let params = &self.params;

        let Xi_key = Uniform::new(-1_f64, 1_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let e = Rq::rand_f64(&mut rng, Xi_err, &params.ring)?;

        let mut s = Rq::rand_f64(&mut rng, Xi_key, &params.ring)?;
        // since s is going to be multiplied by other Rq elements, already
        // compute its NTT
        s.compute_evals();

        let a = Rq::rand_f64(&mut rng, Xi_key, &params.ring)?;

        let pk: PublicKey = PublicKey((&(-a.clone()) * &s) + e, a.clone()); // TODO rm clones
        Ok((SecretKey(s), pk))
    }

    // encrypts a plaintext \in R=Z_Q[X]/(X^N+1)
    fn encrypt(
        &self, // TODO maybe rm?
        mut rng: impl Rng,
        pk: &PublicKey,
        m: &R,
    ) -> Result<(Rq, Rq)> {
        let params = self.params;
        let Xi_key = Uniform::new(-1_f64, 1_f64);
        let Xi_err = Normal::new(0_f64, ERR_SIGMA)?;

        let e_0 = Rq::rand_f64(&mut rng, Xi_err, &params.ring)?;
        let e_1 = Rq::rand_f64(&mut rng, Xi_err, &params.ring)?;

        let v = Rq::rand_f64(&mut rng, Xi_key, &params.ring)?;

        // let m: Rq = Rq::from(*m);
        let m: Rq = m.clone().to_rq(params.ring.q); // TODO rm clone

        Ok((m + e_0 + &v * &pk.0.clone(), &v * &pk.1 + e_1))
    }

    fn decrypt(
        &self, // TODO maybe rm?
        sk: &SecretKey,
        c: (Rq, Rq),
    ) -> Result<R> {
        let m = c.0.clone() + &c.1 * &sk.0;
        Ok(m.mod_centered_q())
    }

    pub fn encode_and_encrypt(
        &self,
        mut rng: impl Rng,
        pk: &PublicKey,
        z: &[C<f64>],
    ) -> Result<(Rq, Rq)> {
        let m: R = self.encoder.encode(&z)?; // polynomial (encoded vec) \in R

        self.encrypt(&mut rng, pk, &m)
    }

    pub fn decrypt_and_decode(&self, sk: SecretKey, c: (Rq, Rq)) -> Result<Vec<C<f64>>> {
        let d = self.decrypt(&sk, c)?;

        self.encoder.decode(&d)
    }

    pub fn add(&self, c0: &(Rq, Rq), c1: &(Rq, Rq)) -> Result<(Rq, Rq)> {
        Ok((&c0.0 + &c1.0, &c0.1 + &c1.1))
    }
    pub fn sub(&self, c0: &(Rq, Rq), c1: &(Rq, Rq)) -> Result<(Rq, Rq)> {
        Ok((&c0.0 - &c1.0, &c0.1 + &c1.1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encrypt_decrypt() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 32;
        let t: u64 = 50;
        let params = Params {
            ring: RingParam { q, n },
            t,
        };
        let scale_factor_u64 = 512_u64; // delta
        let scale_factor = C::<f64>::new(scale_factor_u64 as f64, 0.0); // delta

        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let ckks = CKKS::new(&params, scale_factor);

            let (sk, pk) = ckks.new_key(&mut rng)?;

            let m_raw: R =
                Rq::rand_f64(&mut rng, Uniform::new(0_f64, t as f64), &params.ring)?.to_r();
            let m = &m_raw * &scale_factor_u64;

            let ct = ckks.encrypt(&mut rng, &pk, &m)?;
            let m_decrypted = ckks.decrypt(&sk, ct)?;

            let m_decrypted: Vec<u64> = m_decrypted
                .coeffs()
                .iter()
                .map(|e| (*e as f64 / (scale_factor_u64 as f64)).round() as u64)
                .collect();
            let m_decrypted = Rq::from_vec_u64(&params.ring, m_decrypted);
            // assert_eq!(m_decrypted, Rq::from(m_raw));
            assert_eq!(m_decrypted, m_raw.to_rq(q));
        }

        Ok(())
    }

    #[test]
    fn test_encode_encrypt_decrypt_decode() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 16;
        let t: u64 = 8;
        let params = Params {
            ring: RingParam { q, n },
            t,
        };
        let scale_factor = C::<f64>::new(512.0, 0.0); // delta

        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let ckks = CKKS::new(&params, scale_factor);
            let (sk, pk) = ckks.new_key(&mut rng)?;

            let z: Vec<C<f64>> = std::iter::repeat_with(|| C::<f64>::rand(&mut rng, t))
                .take(n / 2)
                .collect();
            let m: R = ckks.encoder.encode(&z)?;
            println!("{}", m);

            // sanity check
            {
                let z_decoded = ckks.encoder.decode(&m)?;
                let rounded_z_decoded: Vec<C<f64>> = z_decoded
                    .iter()
                    .map(|c| C::<f64>::new(c.re.round(), c.im.round()))
                    .collect();
                assert_eq!(rounded_z_decoded, z);
            }

            let ct = ckks.encrypt(&mut rng, &pk, &m)?;
            let m_decrypted = ckks.decrypt(&sk, ct)?;
            println!("{}", m_decrypted);

            let z_decrypted = ckks.encoder.decode(&m_decrypted)?;

            let rounded_z_decrypted: Vec<C<f64>> = z_decrypted
                .iter()
                .map(|&c| C::<f64>::new(c.re.round(), c.im.round()))
                .collect();
            assert_eq!(rounded_z_decrypted, z);
        }

        Ok(())
    }

    #[test]
    fn test_add() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 16;
        let t: u64 = 8;
        let params = Params {
            ring: RingParam { q, n },
            t,
        };
        let scale_factor = C::<f64>::new(1024.0, 0.0); // delta

        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let ckks = CKKS::new(&params, scale_factor);

            let (sk, pk) = ckks.new_key(&mut rng)?;

            let z0: Vec<C<f64>> = std::iter::repeat_with(|| C::<f64>::rand(&mut rng, t))
                .take(n / 2)
                .collect();
            let z1: Vec<C<f64>> = std::iter::repeat_with(|| C::<f64>::rand(&mut rng, t))
                .take(n / 2)
                .collect();
            let m0: R = ckks.encoder.encode(&z0)?;
            let m1: R = ckks.encoder.encode(&z1)?;

            let ct0 = ckks.encrypt(&mut rng, &pk, &m0)?;
            let ct1 = ckks.encrypt(&mut rng, &pk, &m1)?;

            let ct2 = ckks.add(&ct0, &ct1)?;

            let m2_decrypted = ckks.decrypt(&sk, ct2)?;

            let z_decrypted = ckks.encoder.decode(&m2_decrypted)?;
            let rounded_z_decrypted: Vec<C<f64>> = z_decrypted
                .iter()
                .map(|&c| C::<f64>::new(c.re.round(), c.im.round()))
                .collect();

            let expected_z2: Vec<C<f64>> = itertools::zip_eq(z0, z1).map(|(a, b)| a + b).collect();
            assert_eq!(rounded_z_decrypted, expected_z2);
        }

        Ok(())
    }

    #[test]
    fn test_sub() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 16;
        let t: u64 = 4;
        let params = Params {
            ring: RingParam { q, n },
            t,
        };
        let scale_factor = C::<f64>::new(1024.0, 0.0); // delta

        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let ckks = CKKS::new(&params, scale_factor);

            let (sk, pk) = ckks.new_key(&mut rng)?;

            let z0: Vec<C<f64>> = std::iter::repeat_with(|| C::<f64>::rand(&mut rng, t))
                .take(n / 2)
                .collect();
            let z1: Vec<C<f64>> = std::iter::repeat_with(|| C::<f64>::rand(&mut rng, t))
                .take(n / 2)
                .collect();
            let m0: R = ckks.encoder.encode(&z0)?;
            let m1: R = ckks.encoder.encode(&z1)?;

            let ct0 = ckks.encrypt(&mut rng, &pk, &m0)?;
            let ct1 = ckks.encrypt(&mut rng, &pk, &m1)?;

            let ct2 = ckks.sub(&ct0, &ct1)?;

            let m2_decrypted = ckks.decrypt(&sk, ct2)?;

            let z_decrypted = ckks.encoder.decode(&m2_decrypted)?;
            let rounded_z_decrypted: Vec<C<f64>> = z_decrypted
                .iter()
                .map(|&c| C::<f64>::new(c.re.round(), c.im.round()))
                .collect();

            let expected_z2: Vec<C<f64>> = itertools::zip_eq(z0, z1).map(|(a, b)| a - b).collect();
            assert_eq!(rounded_z_decrypted, expected_z2);
        }

        Ok(())
    }
}
