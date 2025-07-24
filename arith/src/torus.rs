use rand::{distributions::Distribution, Rng};
use std::{
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// Let 𝕋 = ℝ/ℤ, where 𝕋 is a ℤ-module, with homogeneous external product.
/// Let 𝕋q
/// T64 is 𝕋q with q=2^Ω, with Ω=64. We identify 𝕋q=(1/q)ℤ/ℤ ≈ ℤq.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct T64(pub u64);

impl T64 {
    pub fn zero() -> Self {
        Self(0u64)
    }
    pub fn rand(mut rng: impl Rng, dist: impl Distribution<u64>) -> Self {
        let r: u64 = dist.sample(&mut rng);
        Self(r)
    }
    pub fn rand_f64(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        let r: f64 = dist.sample(&mut rng);
        Self(r.round() as u64)
    }

    /// Note: only beta=2 and l=64 is supported.
    pub fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        assert_eq!(beta, 2u32); // only beta=2 supported
        assert_eq!(l, 64u32); // only l=64 supported

        (0..64)
            .rev()
            .map(|i| T64(((self.0 >> i) & 1) as u64))
            .collect()
    }
    pub fn mod_switch<const Q2: u64>(&self) -> T64 {
        todo!()
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
