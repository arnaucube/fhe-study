//! This file implements the struct for an Tuple of Ring Rq elements and its
//! operations, which are performed element-wise.

use anyhow::Result;
use itertools::zip_eq;
use rand::{distributions::Distribution, Rng};
use rand_distr::{Normal, Uniform};
use std::iter::Sum;
use std::{
    array,
    ops::{Add, Mul, Neg, Sub},
};

use crate::Ring;

/// Tuple of K Ring (Rq) elements. We use Vec<R> to allocate it in the heap,
/// since if using a fixed-size array it would overflow the stack.
#[derive(Clone, Debug)]
pub struct TR<R: Ring, const K: usize>(pub Vec<R>);

impl<R: Ring, const K: usize> TR<R, K> {
    pub fn zero() -> Self {
        Self((0..K).into_iter().map(|_| R::zero()).collect())
    }
    pub fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        Self(
            (0..K)
                .into_iter()
                .map(|_| R::rand(&mut rng, &dist))
                .collect(),
        )
    }
    // returns the decomposition of each polynomial element
    pub fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        unimplemented!()
        // self.0.iter().map(|r| r.decompose(beta, l)).collect() // this is Vec<Vec<Vec<R::C>>>
    }
}

impl<R: Ring, const K: usize> TR<R, K> {
    pub fn iter(&self) -> std::slice::Iter<R> {
        self.0.iter()
    }
}

impl<R: Ring, const K: usize> Add<TR<R, K>> for TR<R, K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(
            zip_eq(self.0, other.0)
                .map(|(s, o)| s + o)
                .collect::<Vec<_>>(),
        )
    }
}

impl<R: Ring, const K: usize> Sub<TR<R, K>> for TR<R, K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(zip_eq(self.0, other.0).map(|(s, o)| s - o).collect())
    }
}

impl<R: Ring, const K: usize> Neg for TR<R, K> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.iter().map(|&e_i| -e_i).collect())
    }
}

/// for (TR,TR), the Mul operation is defined as the dot product:
/// for A, B \in R^k, result = Σ A_i * B_i \in R
impl<R: Ring, const K: usize> Mul<TR<R, K>> for TR<R, K> {
    type Output = R;
    fn mul(self, other: Self) -> R {
        zip_eq(self.0, other.0).map(|(s, o)| s * o).sum()
    }
}
impl<R: Ring, const K: usize> Mul<&TR<R, K>> for &TR<R, K> {
    type Output = R;
    fn mul(self, other: &TR<R, K>) -> R {
        zip_eq(self.0.clone(), other.0.clone())
            .map(|(s, o)| s * o)
            .sum()
    }
}

/// for (TR, R), the Mul operation is defined as each element of TR is
/// multiplied by R
impl<R: Ring, const K: usize> Mul<R> for TR<R, K> {
    type Output = TR<R, K>;
    fn mul(self, other: R) -> TR<R, K> {
        Self(self.0.iter().map(|s| s.clone() * other.clone()).collect())
    }
}
impl<R: Ring, const K: usize> Mul<&R> for &TR<R, K> {
    type Output = TR<R, K>;
    fn mul(self, other: &R) -> TR<R, K> {
        TR::<R, K>(self.0.iter().map(|s| s.clone() * other.clone()).collect())
    }
}
