//! This file implements the wrapper on top of the ntt.rs to be able to compute
//! the NTT for non-prime modulus, specifically for modulus 2^64 (for u64).

use crate::ntt::NTT as NTT_p;

// const P0: u64 = 17293822569241362433;

// const P0: u64 = 4611686018427387905;
// const P1: u64 = 4611686018326724609;

// const P0: u64 = 8070449433331580929; // max use 1<<62 -1
// const P1: u64 = 8070450532384645121;

// const P0: u64 = 0x80000000080001; // max use 1<<55 -1
// const P1: u64 = 0x80000000130001;

const BITS: u128 = 64;
const P0: u64 = ((1u128 << BITS) - (1u128 << 28) + 1u128) as u64;
const P1: u64 = ((1u128 << BITS) - (1u128 << 27) + 1u128) as u64;
const P2: u64 = ((1u128 << BITS) - (1u128 << 26) + 1u128) as u64;
// const P0: u64 = ((1u128 << BITS) - (1u128 << 18) + 1u128) as u64;
// const P1: u64 = ((1u128 << BITS) - (1u128 << 17) + 1u128) as u64;
// const P2: u64 = ((1u128 << BITS) - (1u128 << 16) + 1u128) as u64;
// const P0: u64 = 0x0FFFFFFF0000001; // 56 bits each P_i
// const P1: u64 = 0x0FFFFFFE8000001;
// const P2: u64 = 0x0FFFFFFE4000001;

#[derive(Debug)]
pub struct NTT {}

impl NTT {
    pub fn ntt(n: usize, a: &Vec<u64>) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        // apply modulus p_i
        let a_0: Vec<u64> = a.iter().map(|a_i| a_i % P0).collect();
        let a_1: Vec<u64> = a.iter().map(|a_i| a_i % P1).collect();
        let a_2: Vec<u64> = a.iter().map(|a_i| a_i % P2).collect();
        // let a_0: Vec<u64> = a.iter().map(|a_i| (a_i % P0 + P0) % P0).collect();
        // let a_1: Vec<u64> = a.iter().map(|a_i| (a_i % P1 + P1) % P1).collect();
        // let a_2: Vec<u64> = a.iter().map(|a_i| (a_i % P2 + P2) % P2).collect();

        dbg!(&a_0, &a_1, &a_2);
        let r_0 = NTT_p::ntt(P0, n, &a_0);
        let r_1 = NTT_p::ntt(P1, n, &a_1);
        let r_2 = NTT_p::ntt(P2, n, &a_2);
        dbg!(&r_0, &r_1, &r_2);

        (r_0, r_1, r_2)
    }

    pub fn intt(n: usize, r: &(Vec<u64>, Vec<u64>, Vec<u64>)) -> Vec<u64> {
        let a_0 = NTT_p::intt(P0, n, &r.0);
        let a_1 = NTT_p::intt(P1, n, &r.1);
        let a_2 = NTT_p::intt(P2, n, &r.2);
        dbg!(&a_0, &a_1, &a_2);
        // TODO WIP
        // let a_0: Vec<u64> = a_0.iter().map(|a_i| a_i % P0).collect();
        // let a_1: Vec<u64> = a_1.iter().map(|a_i| a_i % P1).collect();
        // let a_2: Vec<u64> = a_2.iter().map(|a_i| a_i % P2).collect();

        reconstruct(a_0, a_1, a_2)
    }
}

/// applies CRT to reconstruct the composite original value
fn reconstruct(a0: Vec<u64>, a1: Vec<u64>, a2: Vec<u64>) -> Vec<u64> {
    let p0: u128 = P0 as u128;
    let p1: u128 = P1 as u128;
    let p2: u128 = P2 as u128;
    // y_i = q/q_i
    let y0: u128 = p1 * p2; // N_i =Q/P0 = P1*P2
    let y1: u128 = p0 * p2;
    let y2: u128 = p0 * p1;

    // y_i^-1 mod q_i = z_i
    dbg!(P0, y0);
    let z0: u128 = inv_mod(p0, y0); // M_i = N_i^-1 mod q_i
    let z1: u128 = inv_mod(p1, y1);
    let z2: u128 = inv_mod(p2, y2);

    // m1 = q1^-1 mod q2
    // aux = (a2 - a1) * m1 mod q2
    // a = a1 + (q1 * m1) * aux

    // m1 == z1
    // let m1 = inv_mod(P1 as u128, P0 as u128); // P0^-1 mod P1
    // let aux: Vec<u64> = itertools::zip_eq(a0.clone(), a1.clone())
    //     .map(|(a0_i, a1_i)| ((a1_i - a0_i) * m1) % P1)
    //     .collect();
    // let a: Vec<u64> = itertools::zip_eq(a0, a1)
    //     // .map(|(a1_i, aux_i)| a1_i + (P1 * m1) * aux_i)
    //     // .map(|(a0_i, aux_i)| a0_i + (P0 * m1) * aux_i)
    //     .map(|(a0_i, a1_i)| a0_i as u128 + ((p0 * z1) % p1) * (((a1_i - a0_i) as u128 * z1) % p1))
    //     .map(|v| v as u64)
    //     // let a: Vec<u64> = itertools::zip_eq(a0, a1)
    //     // .map(|(a0_i, a1_i)| {
    //     //     ((((a0_i as u128 * z0) % P0 as u128) * y0) + (a1_i as u128 * z1 % P1 as u128) * y1)
    //     //         as u64
    //     // })
    //     .collect();
    // a

    // dbg!(a0[0] as u128);
    // dbg!(a0[0] as u128 * y0);
    // dbg!(a0[0] as u128 * z0);
    // dbg!(a0[0] as u128 * y0 * z0);
    dbg!(&y0, &y1, &y2);
    dbg!(&z0, &z1, &z2);
    // WIP, using BigUint to use Q with 192 bits (product of 3 u64)
    let Q = BigUint::from_u64(P0).unwrap()
        * BigUint::from_u64(P1).unwrap()
        * BigUint::from_u64(P2).unwrap();
    let max_u64 = BigUint::from_u128(1_u128 << 64).unwrap();
    // let Q = P0 as u128 * P1 as u128 * P2 as u128;
    let a: Vec<u64> = itertools::multizip((a0, a1, a2))
        .map(|(a0_i, a1_i, a2_i)| {
            // (a0_i as u128 * y0 * z0)// % Q
            //     + (a1_i as u128 * y1 * z1)// % Q
            //     + (a2_i as u128 * y2 * z2) // % Q
            mul_3_big(a0_i as u128, y0, z0, &Q) // TODO rm %Q, since it returns a BigUint, and %Q
                                                // is done later at the end of the additions
                + mul_3_big(a1_i as u128, y1, z1, &Q)
                + mul_3_big(a2_i as u128, y2, z2, &Q)
        })
        .map(|v| {
            // WIP, using BigUint to use Q with 192 bits (product of 3 u64)
            // ((BigUint::from_u128(v).unwrap() % Q.clone()) % max_u64.clone())
            ((v % Q.clone()) % max_u64.clone()).to_u64().unwrap()
        })
        // .map(|v| v % (1 << 63))
        // .map(|v| v % Q)
        // .map(|v| v as u32)
        // .map(|v| v as u64)
        .collect();
    a
    // dbg!(&a);
    // let Q = y2 * P2 as u128;
    // let a: Vec<u128> = a.iter().map(|a_i| a_i % Q).collect();
    // dbg!(&a);
    // let q64 = 1_u128 << 64;
    // let a: Vec<u64> = a.iter().map(|a_i| (a_i % q64) as u64).collect();
    // a
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
/// returns x^-1 mod Q, assuming x and Q are coprime, generally Q is prime
fn inv_mod_new(q: u128, x: u128) -> u128 {
    // by Fermat's Little Theorem, x^-1 mod q \equiv  x^{q-2} mod q
    // exp_mod(q, x, q - 2)
    exp_mod(q, x, q - 2)
}

fn inv_mod(m: u128, a: u128) -> u128 {
    if m == 1 {
        panic!("m==1");
    }

    let mut m = m.clone();
    let mut a = a.clone();
    let m0 = m.clone();
    let mut x0: i128 = 0;
    let mut x1: i128 = 1;

    while a > 1 {
        let q = a / m;
        let t = m.clone();

        m = a % m;
        a = t.clone();

        let t = x0;
        x0 = x1 - (q as i128) * x0;
        x1 = t;
    }

    if x1 < 0 {
        x1 += m0 as i128;
    }

    x1 as u128
}

use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};
fn mul_3_big(a: u128, b: u128, c: u128, Q: &BigUint) -> BigUint {
    let r = (BigUint::from_u128(a).unwrap()
        * BigUint::from_u128(b).unwrap()
        * BigUint::from_u128(c).unwrap());
    // % Q;
    dbg!(&r);
    let r = r % Q;
    dbg!(&r);
    r
    // let max_u64 = BigUint::from_u128(1_u128 << 64).unwrap();
    // dbg!(&r % &max_u64);
    // (r % max_u64).to_u128().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_distr::Distribution;

    use anyhow::Result;

    #[test]
    fn test_inv_mod() -> Result<()> {
        // test vectors used in this test generated in SageMath

        let p: u64 = 0x0FFFFFFF0000001;

        let x = 3;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 48038395846328321);

        let x = (1_u128 << 64) - 1;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 25433122709808565);

        let x = 1_u128 << 120;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 281474976710656);

        // other prime
        let p: u64 = ((1u128 << 64) - (1u128 << 28) + 1u128) as u64;

        let x = 3;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 12297829382294077441);

        let x = (1_u128 << 64) - 1;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 10530692608818076599);

        let x = 1_u128 << 120;
        let x_inv = inv_mod(p as u128, x);
        assert_eq!(x_inv, 18374686616574755070);

        Ok(())
    }

    #[test]
    fn test_reconstruct() -> Result<()> {
        let n: usize = 16;

        use num_bigint::BigUint;
        use num_traits::{FromPrimitive, ToPrimitive};
        let Q = BigUint::from_u64(P0).unwrap()
            * BigUint::from_u64(P1).unwrap()
            * BigUint::from_u64(P2).unwrap();
        let q = BigUint::from_u128(1_u128 << BITS).unwrap();
        let N = BigUint::from_usize(n).unwrap();
        let big2 = BigUint::from_u64(2).unwrap();
        assert!(Q > (N * (&q * &q)) / big2); // sanity check

        use rand::Rng;
        let mut rng = rand::thread_rng();
        let a: Vec<u64> = (0..n)
            .map(|_| rng.gen_range(0..((1u128 << 64) - 1) as u64))
            .collect();
        // let a = vec![14713100818624219214];
        dbg!(&a);

        let a_0: Vec<u64> = a.iter().map(|a_i| a_i % P0).collect();
        let a_1: Vec<u64> = a.iter().map(|a_i| a_i % P1).collect();
        let a_2: Vec<u64> = a.iter().map(|a_i| a_i % P2).collect();
        // let a_0: Vec<u64> = a.iter().map(|a_i| (a_i % P0 + P0) % P0).collect();
        // let a_1: Vec<u64> = a.iter().map(|a_i| (a_i % P1 + P1) % P1).collect();
        // let a_2: Vec<u64> = a.iter().map(|a_i| (a_i % P2 + P2) % P2).collect();
        dbg!(&a_0, &a_1, &a_2);
        let a_reconstructed = reconstruct(a_0, a_1, a_2);

        dbg!(&a_reconstructed);
        assert_eq!(a_reconstructed, a);

        Ok(())
    }

    #[test]
    fn test_ntt() -> Result<()> {
        // println!("{}", 1u128 << 64);
        let n: usize = 1;

        // println!("{}", P0);
        // println!("{}", P1);
        // let q = 1u128 << 64;
        // assert!(P0 as u128 * P1 as u128 > (n as u128 * (q * q)) / 2);

        // let a: Vec<u64> = vec![1u64, 2, 3, 4];
        // let a: Vec<u64> = vec![1u64, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        // let a: Vec<u64> = vec![9u64, 8, 7, 6, 0, 9999, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9];
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let a: Vec<u64> = (0..n)
            .map(|_| rng.gen_range(0..((1u128 << 64) - 1) as u64))
            .collect();

        dbg!(&a);

        let a_ntt = NTT::ntt(n, &a);
        dbg!(&a_ntt);

        let a_intt = NTT::intt(n, &a_ntt);

        dbg!(&a_intt);
        assert_eq!(a_intt, a);

        Ok(())
    }
}
