//! 𝕋_<N,q>[X] = ℝ_<N,q>[X] / ℤ_<N,q>[X], polynomials modulo X^N+1 with
//! coefficients in 𝕋_Q.
//!
//! Note: this is not an algebraic ring, since internal-product is not well
//! defined. But since we work over the discrete torus 𝕋_q, which we identify as
//! 𝕋q = ℤ/qℤ ≈ ℤq, whith q=64. Since we allow product between 𝕋q elements and
//! u64, we fit it into the `Ring` trait (from ring.rs) so that we can compose
//! the 𝕋_<N,q> implementation with the other objects from the code.

use itertools::zip_eq;
use rand::{distributions::Distribution, Rng};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

use crate::{
    ring::{Ring, RingParam},
    torus::T64,
    Rq, Zq,
};

/// 𝕋_<N,Q>[X] = 𝕋<Q>[X]/(X^N +1), polynomials modulo X^N+1 with coefficients in
/// 𝕋, where Q=2^64.
#[derive(Clone, Debug)]
pub struct Tn {
    // pub n: usize,
    pub param: RingParam,
    pub coeffs: Vec<T64>,
}

impl Ring for Tn {
    type C = T64;
    // type Param = usize; // n

    // const Q: u64 = u64::MAX; // WIP
    // const N: usize = N;

    fn param(&self) -> RingParam {
        RingParam {
            q: u64::MAX,
            n: self.param.n,
        }
    }
    fn coeffs(&self) -> Vec<T64> {
        self.coeffs.to_vec()
    }

    fn zero(param: &RingParam) -> Self {
        Self {
            param: *param,
            coeffs: vec![T64::zero(param); param.n],
        }
    }

    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>, param: &RingParam) -> Self {
        Self {
            param: *param,
            coeffs: std::iter::repeat_with(|| T64::rand(&mut rng, &dist, &param))
                .take(param.n)
                .collect(),
        }
        // Self(array::from_fn(|_| T64::rand(&mut rng, &dist)))
    }

    fn from_vec(param: &RingParam, coeffs: Vec<Self::C>) -> Self {
        let mut p = coeffs;
        modulus(param, &mut p);
        Self {
            param: *param,
            coeffs: p,
        }
    }

    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        let elems: Vec<Vec<T64>> = self.coeffs.iter().map(|r| r.decompose(beta, l)).collect();
        // transpose it
        let r: Vec<Vec<T64>> = (0..elems[0].len())
            .map(|i| (0..elems.len()).map(|j| elems[j][i]).collect())
            .collect();
        // convert it to Tn
        r.iter()
            .map(|a_i| Self::from_vec(&self.param, a_i.clone()))
            .collect()
    }

    fn remodule(&self, p: u64) -> Tn {
        todo!()
        // Rq::<P, N>::from_vec_u64(self.coeffs().iter().map(|m_i| m_i.0).collect())
    }

    // fn mod_switch<const P: u64>(&self) -> impl Ring {
    fn mod_switch(&self, p: u64) -> Rq {
        // unimplemented!()
        // TODO WIP
        let coeffs = self
            .coeffs
            .iter()
            .map(|c_i| Zq::from_u64(p, c_i.mod_switch(p).0))
            .collect();
        Rq {
            param: RingParam {
                q: p,
                n: self.param.n,
            },
            coeffs,
            evals: None,
        }
    }

    /// returns [ [(num/den) * self].round() ] mod q
    /// ie. performs the multiplication and division over f64, and then it rounds the
    /// result, only applying the mod Q at the end
    fn mul_div_round(&self, num: u64, den: u64) -> Self {
        let r: Vec<T64> = self
            .coeffs()
            .iter()
            .map(|e| T64(((num as f64 * e.0 as f64) / den as f64).round() as u64))
            .collect();
        Self::from_vec(&self.param, r)
    }
}

impl Tn {
    // multiply self by X^-h
    pub fn left_rotate(&self, h: usize) -> Self {
        let n = self.param.n;

        let h = h % n;
        assert!(h < n);
        let c = &self.coeffs;
        // c[h], c[h+1], c[h+2], ..., c[n-1],  -c[0], -c[1], ..., -c[h-1]
        // let r: Vec<T64> = vec![c[h..N], c[0..h].iter().map(|&c_i| -c_i).collect()].concat();
        let r: Vec<T64> = c[h..n]
            .iter()
            .copied()
            .chain(c[0..h].iter().map(|&x| -x))
            .collect();
        Self::from_vec(&self.param, r)
    }

    pub fn from_vec_u64(param: &RingParam, v: Vec<u64>) -> Self {
        let coeffs = v.iter().map(|c| T64(*c)).collect();
        Self::from_vec(param, coeffs)
    }
}

// apply mod (X^N+1)
pub fn modulus(param: &RingParam, p: &mut Vec<T64>) {
    let n = param.n;
    if p.len() < n {
        return;
    }
    for i in n..p.len() {
        p[i - n] = p[i - n].clone() - p[i].clone();
        p[i] = T64::zero(param);
    }
    p.truncate(n);
}

impl Add<Tn> for Tn {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        // Self(array::from_fn(|i| self.0[i] + rhs.0[i]))
        assert_eq!(self.param, rhs.param);
        Self {
            param: self.param,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l + r)
                .collect(),
        }
    }
}
impl Add<&Tn> for &Tn {
    type Output = Tn;

    fn add(self, rhs: &Tn) -> Self::Output {
        // Tn(array::from_fn(|i| self.0[i] + rhs.0[i]))
        assert_eq!(self.param, rhs.param);
        Tn {
            param: self.param,
            coeffs: zip_eq(self.coeffs.clone(), rhs.coeffs.clone())
                .map(|(l, r)| l + r)
                .collect(),
        }
    }
}
impl AddAssign for Tn {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.param, rhs.param);
        for i in 0..self.param.n {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl Sum<Tn> for Tn {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        // let mut acc = Tn::<N>::zero();
        // for e in iter {
        //     acc += e;
        // }
        // acc
        let first = iter.next().unwrap();
        iter.fold(first, |acc, x| acc + x)
    }
}

impl Sub<Tn> for Tn {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        assert_eq!(self.param, rhs.param);
        Self {
            param: self.param,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l - r)
                .collect(),
        }
    }
}
impl Sub<&Tn> for &Tn {
    type Output = Tn;

    fn sub(self, rhs: &Tn) -> Self::Output {
        // Tn(array::from_fn(|i| self.0[i] - rhs.0[i]))
        assert_eq!(self.param, rhs.param);
        Tn {
            param: self.param,
            coeffs: zip_eq(self.coeffs.clone(), rhs.coeffs.clone())
                .map(|(l, r)| l - r)
                .collect(),
        }
    }
}

impl SubAssign for Tn {
    fn sub_assign(&mut self, rhs: Self) {
        // for i in 0..N {
        //     self.0[i] -= rhs.0[i];
        // }
        assert_eq!(self.param, rhs.param);
        for i in 0..self.param.n {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl Neg for Tn {
    type Output = Self;

    fn neg(self) -> Self::Output {
        // Tn(array::from_fn(|i| -self.0[i]))
        Self {
            param: self.param,
            coeffs: self.coeffs.iter().map(|c_i| -*c_i).collect(),
        }
    }
}

impl PartialEq for Tn {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs && self.param == other.param
    }
}

impl Mul<Tn> for Tn {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        naive_poly_mul(&self, &rhs)
    }
}
impl Mul<&Tn> for &Tn {
    type Output = Tn;

    fn mul(self, rhs: &Tn) -> Self::Output {
        naive_poly_mul(self, rhs)
    }
}

fn naive_poly_mul(poly1: &Tn, poly2: &Tn) -> Tn {
    assert_eq!(poly1.param, poly2.param);
    let n = poly1.param.n;
    let param = poly1.param;

    let poly1: Vec<u128> = poly1.coeffs.iter().map(|c| c.0 as u128).collect();
    let poly2: Vec<u128> = poly2.coeffs.iter().map(|c| c.0 as u128).collect();
    let mut result: Vec<u128> = vec![0; (n * 2) - 1];
    for i in 0..n {
        for j in 0..n {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^n + 1))
    modulus_u128(n, &mut result);

    Tn {
        param,
        // coeffs: array::from_fn(|i| T64(result[i] as u64)),
        coeffs: result.iter().map(|r_i| T64(*r_i as u64)).collect(),
    }
}
fn modulus_u128(n: usize, p: &mut Vec<u128>) {
    if p.len() < n {
        return;
    }
    for i in n..p.len() {
        // p[i - n] = p[i - n].clone() - p[i].clone();
        p[i - n] = p[i - n].wrapping_sub(p[i]);
        p[i] = 0;
    }
    p.truncate(n);
}

impl Mul<T64> for Tn {
    type Output = Self;
    fn mul(self, s: T64) -> Self {
        Self {
            param: self.param,
            // coeffs: array::from_fn(|i| self.coeffs[i] * s),
            coeffs: self.coeffs.iter().map(|c_i| *c_i * s).collect(),
        }
    }
}
// mul by u64
impl Mul<u64> for Tn {
    type Output = Self;
    fn mul(self, s: u64) -> Self {
        // Self(array::from_fn(|i| self.0[i] * s))
        Tn {
            param: self.param,
            coeffs: self.coeffs.iter().map(|c_i| *c_i * s).collect(),
        }
    }
}
impl Mul<&u64> for &Tn {
    type Output = Tn;
    fn mul(self, s: &u64) -> Self::Output {
        // Tn::<N>(array::from_fn(|i| self.0[i] * *s))
        Tn {
            param: self.param,
            coeffs: self.coeffs.iter().map(|c_i| c_i * s).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_left_rotate() {
        let param = RingParam { q: u64::MAX, n: 4 };
        let f = Tn::from_vec(
            &param,
            vec![2i64, 3, -4, -1]
                .iter()
                .map(|c| T64(*c as u64))
                .collect(),
        );

        // expect f*x^-3 == -1 -2x -3x^2 +4x^3
        assert_eq!(
            f.left_rotate(3),
            Tn::from_vec(
                &param,
                vec![-1i64, -2, -3, 4]
                    .iter()
                    .map(|c| T64(*c as u64))
                    .collect(),
            )
        );
        // expect f*x^-1 == 3 -4x -1x^2 -2x^3
        assert_eq!(
            f.left_rotate(1),
            Tn::from_vec(
                &param,
                vec![3i64, -4, -1, -2]
                    .iter()
                    .map(|c| T64(*c as u64))
                    .collect(),
            )
        );
    }
}
