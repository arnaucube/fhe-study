use rand::{distributions::Distribution, Rng};
use std::{
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::ring::Ring;

/// Let 𝕋 = ℝ/ℤ, where 𝕋 is a ℤ-module, with homogeneous external product.
/// Let 𝕋q
/// T64 is 𝕋q with q=2^Ω, with Ω=64. We identify 𝕋q=(1/q)ℤ/ℤ ≈ ℤq.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct T64(pub u64);

// implement the `Ring` trait for T64, so that it can be used where we would use
// `Tn<1>`.
impl Ring for T64 {
    type C = T64;
    const Q: u64 = u64::MAX; // WIP
    const N: usize = 1;

    fn coeffs(&self) -> Vec<T64> {
        vec![self.clone()]
    }
    fn zero() -> Self {
        Self(0u64)
    }
    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        let r: f64 = dist.sample(&mut rng);
        Self(r.round() as u64)
    }
    fn from_vec(coeffs: Vec<Self::C>) -> Self {
        assert_eq!(coeffs.len(), 1);
        coeffs[0]
    }

    // TODO rm beta & l from inputs, make it always beta=2,l=64.
    /// Note: only beta=2 and l=64 is supported.
    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        assert_eq!(beta, 2u32); // only beta=2 supported
                                // assert_eq!(l, 64u32); // only l=64 supported

        // (0..64)
        (0..l)
            .rev()
            .map(|i| T64(((self.0 >> i) & 1) as u64))
            .collect()
    }
    fn remodule<const P: u64>(&self) -> T64 {
        todo!()
    }

    fn mod_switch<const Q2: u64>(&self) -> T64 {
        todo!()
    }

    fn mul_div_round(&self, num: u64, den: u64) -> Self {
        T64(((num as f64 * self.0 as f64) / den as f64).round() as u64)
    }
}

impl T64 {
    pub fn rand_u64(mut rng: impl Rng, dist: impl Distribution<u64>) -> Self {
        let r: u64 = dist.sample(&mut rng);
        Self(r)
    }
}

impl Add<T64> for T64 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0.wrapping_add(rhs.0))
    }
}
impl AddAssign for T64 {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = self.0.wrapping_add(rhs.0)
    }
}

impl Sub<T64> for T64 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0.wrapping_sub(rhs.0))
    }
}
impl SubAssign for T64 {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = self.0.wrapping_sub(rhs.0)
    }
}

impl Neg for T64 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.wrapping_neg())
    }
}

impl Sum for T64 {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self(0), |acc, x| acc + x)
    }
}

impl Mul<T64> for T64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.wrapping_mul(rhs.0))
    }
}
impl MulAssign for T64 {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = self.0.wrapping_mul(rhs.0)
    }
}

// mul by u64
impl Mul<u64> for T64 {
    type Output = Self;

    fn mul(self, s: u64) -> Self {
        if self.0 == 0 {
            return Self(0);
        }
        Self(self.0.wrapping_mul(s))
    }
}
impl Mul<&u64> for &T64 {
    type Output = T64;

    fn mul(self, s: &u64) -> Self::Output {
        T64(self.0.wrapping_mul(*s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::Standard;

    fn recompose(d: Vec<T64>) -> T64 {
        T64(d.iter().fold(0u64, |acc, &b| (acc << 1) | b.0))
    }
    #[test]
    fn test_decompose() {
        let beta: u32 = 2;
        let l: u32 = 64;

        let x = T64(12345);
        let d = x.decompose(beta, l);
        assert_eq!(recompose(d), T64(12345));

        let x = T64(0);
        let d = x.decompose(beta, l);
        assert_eq!(recompose(d), T64(0));

        let x = T64(u64::MAX - 1);
        let d = x.decompose(beta, l);
        assert_eq!(recompose(d), T64(u64::MAX - 1));

        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            let x = T64::rand(&mut rng, Standard);
            let d = x.decompose(beta, l);
            assert_eq!(recompose(d), x);
        }
    }
}
