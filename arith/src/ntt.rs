//! Implementation of the NTT & iNTT, following the CT & GS algorighms, more details in
//! https://eprint.iacr.org/2017/727.pdf, some notes at
//! https://github.com/arnaucube/math/blob/master/notes_ntt.pdf .
//!
//! NOTE: initially I implemented it with fixed Q & N, given as constant
//! generics; but once using real-world parameters, the stack could not handle
//! it, so moved to use Vec instead of fixed-sized arrays, and adapted the NTT
//! implementation to that too.
use crate::{ring::RingParam, ring_nq::Rq, zq::Zq};

use std::collections::HashMap;

#[derive(Debug)]
pub struct NTT {}

use std::sync::{Mutex, OnceLock};

static CACHE: OnceLock<Mutex<HashMap<(u64, usize), (Vec<Zq>, Vec<Zq>, Zq)>>> = OnceLock::new();

fn roots(q: u64, n: usize) -> (Vec<Zq>, Vec<Zq>, Zq) {
    let cache_lock = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    let mut cache = cache_lock.lock().unwrap();
    if let Some(value) = cache.get(&(q, n)) {
        return value.clone();
    }

    let n_inv: Zq = Zq {
        q,
        v: const_inv_mod(q, n as u64),
    };
    let root_of_unity: u64 = primitive_root_of_unity(q, 2 * n);
    let roots_of_unity: Vec<Zq> = roots_of_unity(q, n, root_of_unity);
    let roots_of_unity_inv: Vec<Zq> = roots_of_unity_inv(q, n, roots_of_unity.clone());
    let value = (roots_of_unity, roots_of_unity_inv, n_inv);

    cache.insert((q, n), value.clone());
    value
}

impl NTT {
    /// implements the Cooley-Tukey (CT) algorithm. Details at
    /// https://eprint.iacr.org/2017/727.pdf, also some notes at section 3.1 of
    /// https://github.com/arnaucube/math/blob/master/notes_ntt.pdf
    pub fn ntt(a: &Rq) -> Rq {
        let (q, n) = (a.param.q, a.param.n);
        let (roots_of_unity, _, _) = roots(q, n);

        let mut t = n / 2;
        let mut m = 1;
        let mut r: Vec<Zq> = a.coeffs.clone();
        while m < n {
            let mut k = 0;
            for i in 0..m {
                let S: Zq = roots_of_unity[m + i];
                for j in k..k + t {
                    let U: Zq = r[j];
                    let V: Zq = r[j + t] * S;
                    r[j] = U + V;
                    r[j + t] = U - V;
                }
                k = k + 2 * t;
            }
            t /= 2;
            m *= 2;
        }
        // TODO think if maybe not return a Rq type, or if returned Rq, maybe
        // fill the `evals` field, which is what we're actually returning here
        Rq {
            param: RingParam { q, n },
            coeffs: r,
            evals: None,
        }
    }

    /// implements the Cooley-Tukey (CT) algorithm. Details at
    /// https://eprint.iacr.org/2017/727.pdf, also some notes at section 3.2 of
    /// https://github.com/arnaucube/math/blob/master/notes_ntt.pdf
    pub fn intt(a: &Rq) -> Rq {
        let (q, n) = (a.param.q, a.param.n);
        let (_, roots_of_unity_inv, n_inv) = roots(q, n);

        let mut t = 1;
        let mut m = n / 2;
        let mut r: Vec<Zq> = a.coeffs.clone();
        while m > 0 {
            let mut k = 0;
            for i in 0..m {
                let S: Zq = roots_of_unity_inv[m + i];
                for j in k..k + t {
                    let U: Zq = r[j];
                    let V: Zq = r[j + t];
                    r[j] = U + V;
                    r[j + t] = (U - V) * S;
                }
                k += 2 * t;
            }
            t *= 2;
            m /= 2;
        }
        for i in 0..n {
            r[i] = r[i] * n_inv;
        }
        Rq {
            param: RingParam { q, n },
            coeffs: r,
            // TODO maybe at `evals` place the inputed `a` which is the evals
            // format
            evals: None,
        }
    }
}

/// computes a primitive N-th root of unity using the method described by Thomas
/// Pornin in https://crypto.stackexchange.com/a/63616
const fn primitive_root_of_unity(q: u64, n: usize) -> u64 {
    assert!(n.is_power_of_two());
    assert!((q - 1) % n as u64 == 0);
    let n_u64 = n as u64;

    let mut k = 1;
    while k < q {
        // alternatively could get a random k at each iteration, if so, add the following if:
        // `if k == 0 { continue; }`
        let w = const_exp_mod(q, k, (q - 1) / n_u64);
        if const_exp_mod(q, w, n_u64 / 2) != 1 {
            return w; // w is a primitive N-th root of unity
        }
        k += 1;
    }
    panic!("No primitive root of unity");
}

fn roots_of_unity(q: u64, n: usize, w: u64) -> Vec<Zq> {
    let mut r: Vec<Zq> = vec![Zq { q, v: 0 }; n];
    let mut i = 0;
    let log_n = n.ilog2();
    while i < n {
        // (return the roots in bit-reverset order)
        let j = ((i as u64).reverse_bits() >> (64 - log_n)) as usize;
        r[i] = Zq {
            q,
            v: const_exp_mod(q, w, j as u64),
        };
        i += 1;
    }
    r
}

fn roots_of_unity_inv(q: u64, n: usize, v: Vec<Zq>) -> Vec<Zq> {
    // assumes that the inputted roots are already in bit-reverset order
    let mut r: Vec<Zq> = vec![Zq { q, v: 0 }; n];
    let mut i = 0;
    while i < n {
        r[i] = Zq {
            q,
            v: const_inv_mod(q, v[i].v),
        };
        i += 1;
    }
    r
}

/// returns x^k mod Q
const fn const_exp_mod(q: u64, x: u64, k: u64) -> u64 {
    // work on u128 to avoid overflow
    let mut r = 1u128;
    let mut x = x as u128;
    let mut k = k as u128;
    x = x % q as u128;
    // exponentiation by square strategy
    while k > 0 {
        if k % 2 == 1 {
            r = (r * x) % q as u128;
        }
        x = (x * x) % q as u128;
        k /= 2;
    }
    r as u64
}

/// returns x^-1 mod Q
const fn const_inv_mod(q: u64, x: u64) -> u64 {
    // by Fermat's Little Theorem, x^-1 mod q \equiv  x^{q-2} mod q
    const_exp_mod(q, x, q - 2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Ring;

    use anyhow::Result;

    #[test]
    fn test_ntt() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 4;
        let param = RingParam { q, n };

        let a: Vec<u64> = vec![1u64, 2, 3, 4];
        let a: Rq = Rq::from_vec_u64(&param, a);

        let a_ntt = NTT::ntt(&a);

        let a_intt = NTT::intt(&a_ntt);

        dbg!(&a);
        dbg!(&a_ntt);
        dbg!(&a_intt);
        // dbg!(NTT::ROOT_OF_UNITY);
        // dbg!(NTT::ROOTS_OF_UNITY);

        assert_eq!(a, a_intt);
        Ok(())
    }

    #[test]
    fn test_ntt_loop() -> Result<()> {
        let q: u64 = 2u64.pow(16) + 1;
        let n: usize = 512;
        let param = RingParam { q, n };

        use rand::distributions::Uniform;
        let mut rng = rand::thread_rng();
        let dist = Uniform::new(0_f64, q as f64);

        for _ in 0..1000 {
            let a: Rq = Rq::rand(&mut rng, dist, &param);
            let a_ntt = NTT::ntt(&a);
            let a_intt = NTT::intt(&a_ntt);
            assert_eq!(a, a_intt);
        }
        Ok(())
    }
}
