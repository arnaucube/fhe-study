//! this file implements the non-efficient NTT, which uses multiplication by the
//! Vandermonde matrix.
use crate::zq::Zq;

use anyhow::{anyhow, Result};

#[derive(Debug)]
pub struct NTT<const Q: u64, const N: usize> {
    pub primitive: Zq<Q>,
    // nth_roots: Vec<Zq<Q>>,
    pub ntt: Vec<Vec<Zq<Q>>>,
    pub intt: Vec<Vec<Zq<Q>>>,
}

impl<const Q: u64, const N: usize> NTT<Q, N> {
    pub fn new() -> Result<Self> {
        // TODO change n to be u64 and ensure that is n<Q
        // note: `n` here is not the `N` from `(X^N+1)`
        // TODO: in fact n will be N (trait/struct param)

        // let primitive = Self::get_primitive_root_of_unity((2 * N) as u64)?;
        let primitive = Self::get_primitive_root_of_unity((2 * N) as u64)?;
        // let mut nth_roots = vec![Zq(0); N];
        // let mut w_i = Zq(1);
        // for i in 0..N {
        //     w_i = w_i * primitive;
        //     nth_roots[i] = w_i;
        // }
        let ntt: Vec<Vec<Zq<Q>>> = Self::vandermonde(primitive);
        let intt = Self::invert_vandermonde(&ntt);
        Ok(Self {
            primitive,
            // nth_roots,
            ntt,
            intt,
        })
    }
    /// returns the Vandermonde matrix for the given primitive root of unity.
    /// Vandermonde matrix: https://en.wikipedia.org/wiki/Vandermonde_matrix
    pub fn vandermonde(primitive: Zq<Q>) -> Vec<Vec<Zq<Q>>> {
        let mut v: Vec<Vec<Zq<Q>>> = vec![];
        let n = (2 * N) as u64;
        // let n = N as u64;
        for i in 0..n {
            let mut row: Vec<Zq<Q>> = vec![];
            let primitive_i = primitive.exp(Zq(i));
            let mut primitive_ij = Zq(1);
            for _ in 0..n {
                row.push(primitive_ij);
                primitive_ij = primitive_ij * primitive_i;
            }
            v.push(row);
        }
        v
    }
    // specifically for the Vandermonde matrix
    /// returns the inverse Vandermonde matrix
    pub fn invert_vandermonde(v: &Vec<Vec<Zq<Q>>>) -> Vec<Vec<Zq<Q>>> {
        let n = 2 * N;
        // let n = N;
        let mut inv: Vec<Vec<Zq<Q>>> = vec![];
        for i in 0..n {
            let w_i = v[i][1]; // = w_i^1=w^i^1 = w^i
            let w_i_inv = w_i.inv();
            let mut row: Vec<Zq<Q>> = vec![];
            for j in 0..n {
                row.push(w_i_inv.exp(Zq(j as u64)) / Zq(n as u64));
            }
            inv.push(row);
        }
        inv
    }

    /// computes a primitive N-th root of unity using the method described by
    /// Thomas Pornin in https://crypto.stackexchange.com/a/63616
    pub fn get_primitive_root_of_unity(n: u64) -> Result<Zq<Q>> {
        // using the method described by Thomas Pornin in
        // https://crypto.stackexchange.com/a/63616

        // assert!((Q - 1) % N as u64 == 0);
        assert!((Q - 1) % n == 0);

        // TODO maybe not using Zq and using u64 directly
        let n = Zq(n);
        for k in 0..Q {
            if k == 0 {
                continue;
            }
            let g = Zq(k);
            //  g = F.random_element()
            if g == Zq(0) {
                continue;
            }
            let w = g.exp((-Zq(1)) / n);
            if w.exp(n / Zq(2)) != Zq(1) {
                // g is the generator
                return Ok(w);
            }
        }
        Err(anyhow!("can not find the primitive root of unity"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_distr::Uniform;

    use crate::ringq::matrix_vec_product;
    use crate::ringq::Rq;

    #[test]
    fn roots_of_unity() -> Result<()> {
        const Q: u64 = 12289;
        const N: usize = 512;
        let _ntt = NTT::<Q, N>::new()?;
        Ok(())
    }

    #[test]
    fn vandermonde_ntt() -> Result<()> {
        const Q: u64 = 41;
        const N: usize = 4;
        let primitive = NTT::<Q, N>::get_primitive_root_of_unity((2 * N) as u64)?;
        let v = NTT::<Q, N>::vandermonde(primitive);

        // naively compute the Vandermonde matrix, and assert that the one from the method matches
        // the naively obtained one
        let n2 = (2 * N) as u64;
        let mut v2: Vec<Vec<Zq<Q>>> = vec![];
        for i in 0..n2 {
            let mut row: Vec<Zq<Q>> = vec![];
            for j in 0..n2 {
                row.push(primitive.exp(Zq(i * j)));
            }
            v2.push(row);
        }
        assert_eq!(v, v2);

        let v_inv = NTT::<Q, N>::invert_vandermonde(&v);

        let mut rng = rand::thread_rng();
        let uniform_distr = Uniform::new(0_f64, Q as f64);
        let a = Rq::<Q, N>::rand_f64(&mut rng, uniform_distr)?;
        // let a = PR::<Q, N>::new_from_u64(vec![36, 21, 9, 19]);

        // let a_padded_coeffs: [Zq<Q>; 2 * N] =
        //     std::array::from_fn(|i| if i < N { a.coeffs[i] } else { Zq::zero() });
        let mut a_padded = a.coeffs.to_vec();
        a_padded.append(&mut vec![Zq(0); N]);
        // let a_ntt = a_padded.mul_by_matrix(&v)?;
        let a_ntt = matrix_vec_product(&v, &a_padded)?;
        let a_intt: Vec<Zq<Q>> = matrix_vec_product(&v_inv, &a_ntt)?;
        assert_eq!(a_intt, a_padded);
        let a_intt_arr: [Zq<Q>; N] = std::array::from_fn(|i| a_intt[i]);
        assert_eq!(Rq::new(a_intt_arr, None), a);

        Ok(())
    }

    #[test]
    fn vec_by_ntt() -> Result<()> {
        const Q: u64 = 257;
        const N: usize = 4;
        // let primitive = NTT::<Q, N>::get_primitive_root_of_unity((2*N) as u64)?;
        let ntt = NTT::<Q, N>::new()?;

        let a: Vec<Zq<Q>> = vec![256, 256, 256, 256, 0, 0, 0, 0]
            .iter()
            .map(|&e| Zq::from_u64(e))
            .collect();
        let a_ntt = matrix_vec_product(&ntt.ntt, &a)?;
        let a_intt = matrix_vec_product(&ntt.intt, &a_ntt)?;
        assert_eq!(a_intt, a);

        Ok(())
    }

    #[test]
    fn bench_ntt() -> Result<()> {
        // const Q: u64 = 12289;
        // const N: usize = 512;
        const Q: u64 = 257;
        const N: usize = 4;
        // let primitive = NTT::<Q, N>::get_primitive_root_of_unity((2*N) as u64)?;
        let ntt = NTT::<Q, N>::new()?;

        let rng = rand::thread_rng();
        let a = Rq::<Q, { 2 * N }>::rand_f64(rng, Uniform::new(0_f64, (Q - 1) as f64))?;
        let a = a.coeffs;
        dbg!(&a);
        let a_ntt = matrix_vec_product(&ntt.ntt, &a.to_vec())?;
        dbg!(&a_ntt);
        let a_intt = matrix_vec_product(&ntt.intt, &a_ntt)?;
        dbg!(&a_intt);
        assert_eq!(a_intt, a);
        // TODO bench

        Ok(())
    }
}
