//! Implementation of the NTT & iNTT, following the CT & GS algorighms, more
//! details in https://github.com/arnaucube/math/blob/master/notes_ntt.pdf .
use crate::zq::Zq;

#[derive(Debug)]
pub struct NTT<const Q: u64, const N: usize> {}

impl<const Q: u64, const N: usize> NTT<Q, N> {
    const N_INV: Zq<Q> = Zq(const_inv_mod::<Q>(N as u64));
    // since we work over Zq[X]/(X^N+1) (negacyclic), get the 2*N-th root of unity
    pub(crate) const ROOT_OF_UNITY: u64 = primitive_root_of_unity::<Q>(2 * N);
    pub(crate) const ROOTS_OF_UNITY: [Zq<Q>; N] = roots_of_unity(Self::ROOT_OF_UNITY);
    const ROOTS_OF_UNITY_INV: [Zq<Q>; N] = roots_of_unity_inv(Self::ROOTS_OF_UNITY);
}

impl<const Q: u64, const N: usize> NTT<Q, N> {
    /// implements the Cooley-Tukey (CT) algorithm. Details at section 3.1 of
    /// https://github.com/arnaucube/math/blob/master/notes_ntt.pdf
    pub fn ntt(a: [Zq<Q>; N]) -> [Zq<Q>; N] {
        let mut t = N / 2;
        let mut m = 1;
        let mut r: [Zq<Q>; N] = a.clone();
        while m < N {
            let mut k = 0;
            for i in 0..m {
                let S: Zq<Q> = Self::ROOTS_OF_UNITY[m + i];
                for j in k..k + t {
                    let U: Zq<Q> = r[j];
                    let V: Zq<Q> = r[j + t] * S;
                    r[j] = U + V;
                    r[j + t] = U - V;
                }
                k = k + 2 * t;
            }
            t /= 2;
            m *= 2;
        }
        r
    }

    /// implements the Gentleman-Sande (GS) algorithm. Details at section 3.2 of
    /// https://github.com/arnaucube/math/blob/master/notes_ntt.pdf
    pub fn intt(a: [Zq<Q>; N]) -> [Zq<Q>; N] {
        let mut t = 1;
        let mut m = N / 2;
        let mut r: [Zq<Q>; N] = a.clone();
        while m > 0 {
            let mut k = 0;
            for i in 0..m {
                let S: Zq<Q> = Self::ROOTS_OF_UNITY_INV[m + i];
                for j in k..k + t {
                    let U: Zq<Q> = r[j];
                    let V: Zq<Q> = r[j + t];
                    r[j] = U + V;
                    r[j + t] = (U - V) * S;
                }
                k += 2 * t;
            }
            t *= 2;
            m /= 2;
        }
        for i in 0..N {
            r[i] = r[i] * Self::N_INV;
        }
        r
    }
}

/// computes a primitive N-th root of unity using the method described by Thomas
/// Pornin in https://crypto.stackexchange.com/a/63616
const fn primitive_root_of_unity<const Q: u64>(N: usize) -> u64 {
    assert!(N.is_power_of_two());
    assert!((Q - 1) % N as u64 == 0);

    let n: u64 = N as u64;
    let mut k = 1;
    while k < Q {
        // alternatively could get a random k at each iteration, if so, add the following if:
        // `if k == 0 { continue; }`
        let w = const_exp_mod::<Q>(k, (Q - 1) / n);
        if const_exp_mod::<Q>(w, n / 2) != 1 {
            return w; // w is a primitive N-th root of unity
        }
        k += 1;
    }
    panic!("No primitive root of unity");
}

const fn roots_of_unity<const Q: u64, const N: usize>(w: u64) -> [Zq<Q>; N] {
    let mut r: [Zq<Q>; N] = [Zq(0u64); N];
    let mut i = 0;
    let log_n = N.ilog2();
    while i < N {
        // (return the roots in bit-reverset order)
        let j = ((i as u64).reverse_bits() >> (64 - log_n)) as usize;
        r[i] = Zq(const_exp_mod::<Q>(w, j as u64));
        i += 1;
    }
    r
}

const fn roots_of_unity_inv<const Q: u64, const N: usize>(v: [Zq<Q>; N]) -> [Zq<Q>; N] {
    // assumes that the inputted roots are already in bit-reverset order
    let mut r: [Zq<Q>; N] = [Zq(0u64); N];
    let mut i = 0;
    while i < N {
        r[i] = Zq(const_inv_mod::<Q>(v[i].0));
        i += 1;
    }
    r
}

/// returns x^k mod Q
const fn const_exp_mod<const Q: u64>(x: u64, k: u64) -> u64 {
    let mut r = 1u64;
    let mut x = x;
    let mut k = k;
    x = x % Q;
    // exponentiation by square strategy
    while k > 0 {
        if k % 2 == 1 {
            r = (r * x) % Q;
        }
        x = (x * x) % Q;
        k /= 2;
    }
    r
}

/// returns x^-1 mod Q
const fn const_inv_mod<const Q: u64>(x: u64) -> u64 {
    // by Fermat's Little Theorem, x^-1 mod q \equiv  x^{q-2} mod q
    const_exp_mod::<Q>(x, Q - 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    use anyhow::Result;
    use std::array;

    #[test]
    fn test_ntt() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 4;

        let a: [u64; N] = [1u64, 2, 3, 4];
        let a: [Zq<Q>; N] = array::from_fn(|i| Zq::new(a[i]));

        let a_ntt = NTT::<Q, N>::ntt(a);

        let a_intt = NTT::<Q, N>::intt(a_ntt);

        dbg!(&a);
        dbg!(&a_ntt);
        dbg!(&a_intt);
        dbg!(NTT::<Q, N>::ROOT_OF_UNITY);
        dbg!(NTT::<Q, N>::ROOTS_OF_UNITY);

        assert_eq!(a, a_intt);
        Ok(())
    }

    #[test]
    fn test_ntt_loop() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 512;

        use rand::distributions::Distribution;
        use rand::distributions::Uniform;
        let mut rng = rand::thread_rng();
        let dist = Uniform::new(0_f64, Q as f64);

        for _ in 0..100 {
            let a: [Zq<Q>; N] = array::from_fn(|_| Zq::from_f64(dist.sample(&mut rng)));
            let a_ntt = NTT::<Q, N>::ntt(a);
            let a_intt = NTT::<Q, N>::intt(a_ntt);
            assert_eq!(a, a_intt);
        }
        Ok(())
    }
}
