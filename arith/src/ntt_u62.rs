//! This file implements the wrapper on top of the ntt.rs to be able to compute
//! the NTT for non-prime modulus, specifically for modulus 2^64 (for u64).

use crate::ntt::NTT as NTT_p;

// const P0: u64 = 17293822569241362433;

// const P0: u64 = 4611686018427387905;
// const P1: u64 = 4611686018326724609;

const P0: u64 = 8070449433331580929; // max use 1<<62
const P1: u64 = 8070450532384645121;

// const P0: u64 = 0x80000000080001; // max use 1<<55
// const P1: u64 = 0x80000000130001;

// const P0: u64 = ((1u128 << 64) - (1u128 << 28) + 1u128) as u64;
// const P1: u64 = ((1u128 << 64) - (1u128 << 27) + 1u128) as u64;
// const P0: u64 = ((1u128 << 64) - (1u128 << 28) + 1u128) as u64;
// const P1: u64 = ((1u128 << 64) - (1u128 << 27) + 1u128) as u64;
// const P2: u64 = (1 << 60) - (1 << 26) + 1;

#[derive(Debug)]
pub struct NTT {}

impl NTT {
    pub fn ntt(
        n: usize,
        a: &Vec<u64>,
    ) -> (
        Vec<u64>,
        Vec<u64>,
        // Vec<u64>,
        // Vec<u64>,
        // Vec<u64>,
        // Vec<u64>,
        // Vec<u64>,
    ) {
        // TODO ensure that: a_i <P0

        // apply modulus p_i
        let a_0: Vec<u64> = a.iter().map(|a_i| (a_i % P0 + P0) % P0).collect();
        let a_1: Vec<u64> = a.iter().map(|a_i| (a_i % P1 + P1) % P1).collect();
        // let a_2: Vec<u64> = a.iter().map(|a_i| (a_i % P2 + P2) % P2).collect();
        // let a_3: Vec<u64> = a.iter().map(|a_i| (a_i % P3 + P3) % P3).collect();
        // let a_4: Vec<u64> = a.iter().map(|a_i| (a_i % P4 + P4) % P4).collect();
        // let a_5: Vec<u64> = a.iter().map(|a_i| (a_i % P5 + P5) % P5).collect();
        // let a_6: Vec<u64> = a.iter().map(|a_i| (a_i % P6 + P6) % P6).collect();

        let r_0 = NTT_p::ntt(P0, n, &a_0);
        let r_1 = NTT_p::ntt(P1, n, &a_1);
        // let r_2 = NTT_p::ntt(P2, n, &a_2);
        // let r_3 = NTT_p::ntt(P3, n, &a_3);
        // let r_4 = NTT_p::ntt(P4, n, &a_4);
        // let r_5 = NTT_p::ntt(P5, n, &a_5);
        // let r_6 = NTT_p::ntt(P6, n, &a_6);

        (r_0, r_1) //, r_2) //, r_3, r_4) //, r_5, r_6)
    }

    pub fn intt(
        n: usize,
        r: &(
            Vec<u64>,
            Vec<u64>,
            // Vec<u64>,
            // Vec<u64>,
            // Vec<u64>,
            // Vec<u64>,
            // Vec<u64>,
        ),
    ) -> Vec<u64> {
        let a_0 = NTT_p::intt(P0, n, &r.0);
        let a_1 = NTT_p::intt(P1, n, &r.1);
        // let a_2 = NTT_p::intt(P2, n, &r.2);
        // let a_3 = NTT_p::intt(P3, n, &r.3);
        // let a_4 = NTT_p::intt(P4, n, &r.4);
        // let a_5 = NTT_p::intt(P5, n, &r.5);
        // let a_6 = NTT_p::intt(P6, n, &r.6);

        // Garner CRT for two moduli: combine (r1 mod p1, r2 mod p2) -> Z/(p1*p2)
        // let inv_p1_mod_p2: u128 = inv_mod_u64(p1 % p2, p2) as u128;
        // const INV_P1_MOD_P2: u128 = 4895217125691974194;

        reconstruct(a_0, a_1) //, a_2) // , a_3, a_4) //, a_5, a_6)
    }
}

fn reconstruct(
    a0: Vec<u64>,
    a1: Vec<u64>,
    // a2: Vec<u64>,
    // a_3: Vec<u64>,
    // a_4: Vec<u64>,
    // a_5: Vec<u64>,
    // a_6: Vec<u64>,
) -> Vec<u64> {
    // let Q = P0 as u128 * P1 as u128;

    // y_i = q/q_i
    // let y0 = ((u64::MAX as u128 + 1) / P0 as u128);
    // let y1 = ((u64::MAX as u128 + 1) / P1 as u128);
    // let y2 = ((u64::MAX as u128 + 1) / P2 as u128) as u64;
    let y0: u128 = P1 as u128; // N_i =Q/P0 = P1*P2
    let y1: u128 = P0 as u128;
    // let y1: u128 = P0 as u128 * P2 as u128;
    // let y0: u128 = P1 as u128 * P2 as u128; // N_i =Q/P0 = P1*P2
    // let y1: u128 = P0 as u128 * P2 as u128;
    // let y2: u128 = P0 as u128 * P1 as u128;
    // let y0 = (Q / P0 as u128) as u64;
    // let y1 = (Q / P1 as u128) as u64;
    // let y2 = ((u64::MAX as u128 + 1) / P2 as u128) as u64;
    // let y3 = ((u64::MAX as u128 + 1) / P3 as u128) as u64;
    // let y4 = ((u64::MAX as u128 + 1) / P4 as u128) as u64;

    // y_i^-1 mod q_i = z_i
    let z0: u128 = inv_mod(P0 as u128, y0); // M_i = N_i^-1 mod q_i
    let z1: u128 = inv_mod(P1 as u128, y1);
    // let z2: u128 = inv_mod(P2 as u128, y2);
    // let y2_inv = inv_mod(P2 as u128, y2);
    // let y3_inv = inv_mod(P3 as u128, y3);
    // let y4_inv = inv_mod(P4 as u128, y4);

    // m1 = q1^-1 mod q2
    // aux = (a2 - a1) * m1 mod q2
    // a = a1 + (q1 * m1) * aux

    /*
    let m1 = inv_mod(P1 as u128, P0 as u128) as u64; // P0^-1 mod P1
    let aux: Vec<u64> = itertools::zip_eq(a0.clone(), a1.clone())
        .map(|(a0_i, a1_i)| ((a1_i - a0_i) * m1) % P1)
        .collect();
    let a: Vec<u64> = itertools::zip_eq(a0, aux)
        // .map(|(a1_i, aux_i)| a1_i + (P1 * m1) * aux_i)
        // .map(|(a0_i, aux_i)| a0_i + (P0 * m1) * aux_i)
        .map(|(a0_i, aux_i)| a0_i + ((P0 * m1) % P1) * aux_i)
        .collect();
    a
        */
    let p0: u128 = P0 as u128;
    let p1: u128 = P1 as u128;
    let a: Vec<u64> = itertools::zip_eq(a0, a1)
        .map(|(a0_i, a1_i)| a0_i as u128 + ((p0 * z1) % p1) * (((a1_i - a0_i) as u128 * z1) % p1))
        .map(|v| v as u64)
        .collect();
    a
    // dbg!(a0[0] as u128);
    // dbg!(a0[0] as u128 * y0);
    // dbg!(a0[0] as u128 * z0);
    // dbg!(a0[0] as u128 * y0 * z0);
    // let a: Vec<u128> = itertools::multizip((a0, a1, a2))
    //     .map(|(a0_i, a1_i, a2_i)| {
    //         a0_i as u128 * y0 * z0 + a1_i as u128 * y1 * z1 + a2_i as u128 * y2 * z2
    //     })
    //     .collect();
    // dbg!(&a);
    // let Q = y2 * P2 as u128;
    // let a: Vec<u128> = a.iter().map(|a_i| a_i % Q).collect();
    // dbg!(&a);
    // let q64 = 1_u128 << 64;
    // let a: Vec<u64> = a.iter().map(|a_i| (a_i % q64) as u64).collect();
    // a

    /*
    // x_i*z_i mod q_i
    let r0: Vec<u64> = a_0.iter().map(|a_i| ((a_i * z0) % P0) * y0).collect();
    let r1: Vec<u64> = a_1.iter().map(|a_i| ((a_i * z1) % P1) * y1).collect();
    // let r0: Vec<u64> = a_0.iter().map(|a_i| ((a_i * z0) % P0) * y0).collect();
    // let r1: Vec<u64> = a_1.iter().map(|a_i| ((a_i * z1) % P1) * y1).collect();
    // let r2: Vec<u64> = a_2.iter().map(|a_i| ((a_i * y2_inv) % P2) * y2).collect();
    // let r3: Vec<u64> = a_3.iter().map(|a_i| ((a_i * y3_inv) % P3) * y3).collect();
    // let r4: Vec<u64> = a_4.iter().map(|a_i| ((a_i * y4_inv) % P4) * y4).collect();

    let r: Vec<u64> = itertools::multizip((r0.iter(), r1.iter()))
        .map(|(a, b)| a + b)
        .collect();
    // let r = r0;
    //
    dbg!(&r);

    let p1p2: u128 = (P0 as u128) * (P1 as u128);
    // let p1p2_inv: u128 = inv_mod((P0 % P1) as u128, P1) as u128;
    let p1p2_inv: u128 = inv_mod((P0) as u128, P1) as u128;
    dbg!(&p1p2);
    dbg!(&p1p2_inv);
    // let p1p2: u128 = P0 as u128 / 2; // PIHALF
    let r = r
        .iter()
        .map(|c_i_u64| {
            let c_i = *c_i_u64 as u128;
            if c_i * 2 >= p1p2 {
                // if c_i >= p1p2 {
                c_i.wrapping_sub(p1p2) as u64
            } else {
                c_i as u64
            }
        })
        .collect();
    // let r: Vec<u64> = itertools::multizip((r0.iter(), r1.iter(), r2.iter(), r3.iter(), r4.iter()))
    //     .map(|(a, b, c, d, e)| a + b + c + d + e)
    //     .collect();
    // let mut r = a_0 + y0_inv + a_1 * y1_inv + a_2 * y2_inv + a_3 * y3_inv + a_4 * y4_inv;

    r
    */
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
    use rand_distr::Distribution;

    use anyhow::Result;

    #[test]
    fn test_dbg() -> Result<()> {
        println!("{}", 1u128 << 64);
        let n: usize = 16;

        println!("{}", P0);
        println!("{}", P1);
        // let q = 1u128 << 64;
        // assert!(P0 as u128 * P1 as u128 > (n as u128 * (q * q)) / 2);

        // let a: Vec<u64> = vec![1u64, 2, 3, 4];
        // let a: Vec<u64> = vec![1u64, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        // let a: Vec<u64> = vec![9u64, 8, 7, 6, 0, 9999, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9];
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let a: Vec<u64> = (0..n)
            // .map(|_| rng.gen_range(0..=(1u64 << 57) - 1) - (1u64 << 56))
            .map(|_| rng.gen_range(0..(1 << 62)))
            .collect();

        dbg!(a.len());

        let a_ntt = NTT::ntt(n, &a);
        dbg!(&a_ntt);

        let a_intt = NTT::intt(n, &a_ntt);

        dbg!(&a_intt);
        assert_eq!(a_intt, a);

        Ok(())
    }
}
