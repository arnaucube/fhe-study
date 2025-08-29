//! This file implements the wrapper on top of the ntt.rs to be able to compute
//! the NTT for non-prime modulus, specifically for modulus 2^64 (for u64).

use crate::ntt::NTT as NTT_p;

const P0: u64 = 8070449433331580929; // max use: 1<<62
const P1: u64 = 8070450532384645121;

// const P0: u64 = 0x80000000080001; // max use: 1<<55
// const P1: u64 = 0x80000000130001;

#[derive(Debug)]
pub struct NTT {}

impl NTT {
    pub fn ntt(n: usize, a: &Vec<u64>) -> (Vec<u64>, Vec<u64>) {
        // apply modulus p_i
        let a_0: Vec<u64> = a.iter().map(|a_i| (a_i % P0 + P0) % P0).collect();
        let a_1: Vec<u64> = a.iter().map(|a_i| (a_i % P1 + P1) % P1).collect();

        let r_0 = NTT_p::ntt(P0, n, &a_0);
        let r_1 = NTT_p::ntt(P1, n, &a_1);

        (r_0, r_1)
    }

    pub fn intt(n: usize, r: &(Vec<u64>, Vec<u64>)) -> Vec<u64> {
        let a_0 = NTT_p::intt(P0, n, &r.0);
        let a_1 = NTT_p::intt(P1, n, &r.1);

        reconstruct(a_0, a_1)
    }
}

/// applies CRT to reconstruct the composite original value. Uses Garner's CRT algorithm for two
/// moduli: combine (r1 mod p1, r2 mod p2) -> Z/(p1*p2)
fn reconstruct(a0: Vec<u64>, a1: Vec<u64>) -> Vec<u64> {
    // y_i = q/q_i
    // let y0: u128 = P1 as u128;
    let y1: u128 = P0 as u128; // Q/P1 = P0*P1/P1 = P0

    // y_i^-1 mod q_i = z_i
    // let z0: u128 = inv_mod(P0 as u128, y0); // M_i = N_i^-1 mod q_i
    let z1: u128 = inv_mod(P1 as u128, y1);

    // m1 = q1^-1 mod q2
    // aux = (a2 - a1) * m1 mod q2
    // a = a1 + (q1 * m1) * aux
    let p0: u128 = P0 as u128;
    let p1: u128 = P1 as u128;
    let a: Vec<u64> = itertools::zip_eq(a0, a1)
        .map(|(a0_i, a1_i)| a0_i as u128 + ((p0 * z1) % p1) * (((a1_i - a0_i) as u128 * z1) % p1))
        .map(|v| v as u64)
        .collect();
    a
}

fn exp_mod(q: u128, x: u128, k: u128) -> u128 {
    // work on u128 to avoid overflow
    let mut r = 1u128;
    let mut x = x.clone();
    let mut k = k.clone();
    x = x % q;
    // exponentiation by square strategy
    while k > 0 {
        if k % 2 == 1 {
            r = (r * x) % q;
        }
        x = (x * x) % q;
        k /= 2;
    }
    r
}
/// returns x^-1 mod Q
fn inv_mod(q: u128, x: u128) -> u128 {
    // by Fermat's Little Theorem, x^-1 mod q \equiv  x^{q-2} mod q
    exp_mod(q, x, q - 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    use anyhow::Result;

    #[test]
    fn test_dbg() -> Result<()> {
        let n: usize = 16;

        // let q = 1u128 << 64;
        // assert!(P0 as u128 * P1 as u128 > (n as u128 * (q * q)) / 2);

        use rand::Rng;
        let mut rng = rand::thread_rng();
        let a: Vec<u64> = (0..n).map(|_| rng.gen_range(0..(1 << 62))).collect();

        dbg!(a.len());

        let a_ntt = NTT::ntt(n, &a);
        dbg!(&a_ntt);

        let a_intt = NTT::intt(n, &a_ntt);

        dbg!(&a_intt);
        assert_eq!(a_intt, a);

        Ok(())
    }
}
