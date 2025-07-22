//! 𝕋_<N,q>[X] = ℝ_<N,q>[X] / ℤ_<N,q>[X], polynomials modulo X^N+1 with
//! coefficients in 𝕋_Q.
//!
//! Note: this is not an algebraic ring, since internal-product is not well
//! defined. But since we work over the discrete torus 𝕋_q, which we identify as
//! 𝕋q = ℤ/qℤ ≈ ℤq, whith q=64. Since we allow product between 𝕋q elements and
//! u64, we fit it into the `Ring` trait (from ring.rs) so that we can compose
//! the 𝕋_<N,q> implementation with the other objects from the code.

use rand::{distributions::Distribution, Rng};
use std::array;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

use crate::{ring::Ring, torus::T64};

/// 𝕋_<N,Q>[X] = 𝕋<Q>[X]/(X^N +1), polynomials modulo X^N+1 with coefficients in
/// 𝕋, where Q=2^64.
#[derive(Clone, Copy, Debug)]
pub struct Tn<const N: usize>(pub [T64; N]);

impl<const N: usize> Ring for Tn<N> {
    type C = T64;

    fn coeffs(&self) -> Vec<T64> {
        self.0.to_vec()
    }

    fn zero() -> Self {
        Self(array::from_fn(|_| T64::zero()))
    }

    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        Self(array::from_fn(|_| T64::rand_f64(&mut rng, &dist)))
    }

    fn from_vec(coeffs: Vec<Self::C>) -> Self {
        let mut p = coeffs;
        modulus::<N>(&mut p);
        Self(array::from_fn(|i| p[i]))
    }

    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        let elems: Vec<Vec<T64>> = self.0.iter().map(|r| r.decompose(beta, l)).collect();
        // transpose it
        let r: Vec<Vec<T64>> = (0..elems[0].len())
            .map(|i| (0..elems.len()).map(|j| elems[j][i]).collect())
            .collect();
        // convert it to Tn<N>
        r.iter().map(|a_i| Self::from_vec(a_i.clone())).collect()
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
        Self::from_vec(r)
    }
}

// apply mod (X^N+1)
pub fn modulus<const N: usize>(p: &mut Vec<T64>) {
    if p.len() < N {
        return;
    }
    for i in N..p.len() {
        p[i - N] = p[i - N].clone() - p[i].clone();
        p[i] = T64::zero();
    }
    p.truncate(N);
}

impl<const N: usize> Add<Tn<N>> for Tn<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}
impl<const N: usize> Add<&Tn<N>> for &Tn<N> {
    type Output = Tn<N>;

    fn add(self, rhs: &Tn<N>) -> Self::Output {
        Tn(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}
impl<const N: usize> AddAssign for Tn<N> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.0[i] += rhs.0[i];
        }
    }
}

impl<const N: usize> Sum<Tn<N>> for Tn<N> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = Tn::<N>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<const N: usize> Sub<Tn<N>> for Tn<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}
impl<const N: usize> Sub<&Tn<N>> for &Tn<N> {
    type Output = Tn<N>;

    fn sub(self, rhs: &Tn<N>) -> Self::Output {
        Tn(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}

impl<const N: usize> SubAssign for Tn<N> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.0[i] -= rhs.0[i];
        }
    }
}

impl<const N: usize> PartialEq for Tn<N> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<const N: usize> Mul<Tn<N>> for Tn<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        naive_poly_mul(&self, &rhs)
    }
}
impl<const N: usize> Mul<&Tn<N>> for &Tn<N> {
    type Output = Tn<N>;

    fn mul(self, rhs: &Tn<N>) -> Self::Output {
        naive_poly_mul(self, rhs)
    }
}

fn naive_poly_mul<const N: usize>(poly1: &Tn<N>, poly2: &Tn<N>) -> Tn<N> {
    let poly1: Vec<u128> = poly1.0.iter().map(|c| c.0 as u128).collect();
    let poly2: Vec<u128> = poly2.0.iter().map(|c| c.0 as u128).collect();
    let mut result: Vec<u128> = vec![0; (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    modulus_u128::<N>(&mut result);

    // sanity check: check that there are no coeffs > i64_max
    assert_eq!(
        result,
        Tn::<N>(array::from_fn(|i| T64(result[i] as u64)))
            .coeffs()
            .iter()
            .map(|c| c.0 as u128)
            .collect::<Vec<_>>()
    );
    Tn(array::from_fn(|i| T64(result[i] as u64)))
}
fn modulus_u128<const N: usize>(p: &mut Vec<u128>) {
    if p.len() < N {
        return;
    }
    for i in N..p.len() {
        p[i - N] = p[i - N].clone() - p[i].clone();
        p[i] = 0;
    }
    p.truncate(N);
}

// mul by u64
impl<const N: usize> Mul<u64> for Tn<N> {
    type Output = Self;
    fn mul(self, s: u64) -> Self {
        Self(array::from_fn(|i| self.0[i] * s))
    }
}
impl<const N: usize> Mul<&u64> for &Tn<N> {
    type Output = Tn<N>;
    fn mul(self, s: &u64) -> Self::Output {
        Tn::<N>(array::from_fn(|i| self.0[i] * *s))
    }
}
