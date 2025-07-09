//! This file implements the struct for an Tuple of Ring Rq elements and its
//! operations.

use anyhow::Result;
use itertools::zip_eq;
use rand::{distributions::Distribution, Rng};
use rand_distr::{Normal, Uniform};
use std::iter::Sum;
use std::{array, ops};

use crate::Ring;

// #[derive(Clone, Copy, Debug)]
// pub struct TR<R: Ring, const K: usize>([R; K]);

/// Tuple of K Ring (Rq) elements. We use Vec<R> to allocate it in the heap,
/// since if using a fixed-size array it would overflow the stack.
#[derive(Clone, Debug)]
pub struct TR<R: Ring, const K: usize>(Vec<R>);

impl<R: Ring, const K: usize> TR<R, K> {
    pub fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        Self(
            (0..K)
                .into_iter()
                .map(|_| R::rand(&mut rng, &dist))
                .collect(),
        )
    }
}

impl<R: Ring, const K: usize> ops::Add<TR<R, K>> for TR<R, K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(
            zip_eq(self.0, other.0)
                .map(|(s, o)| s + o)
                .collect::<Vec<_>>(),
        )
    }
}

impl<R: Ring, const K: usize> ops::Sub<TR<R, K>> for TR<R, K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(zip_eq(self.0, other.0).map(|(s, o)| s - o).collect())
    }
}

/// for (TR,TR), the Mul operation is defined as:
/// for A, B \in R^k, result = Σ A_i * B_i \in R
impl<R: Ring, const K: usize> ops::Mul<TR<R, K>> for TR<R, K> {
    type Output = R;
    fn mul(self, other: Self) -> R {
        zip_eq(self.0, other.0).map(|(s, o)| s * o).sum()
    }
}
impl<R: Ring, const K: usize> ops::Mul<&TR<R, K>> for &TR<R, K> {
    type Output = R;
    fn mul(self, other: &TR<R, K>) -> R {
        zip_eq(self.0.clone(), other.0.clone())
            .map(|(s, o)| s * o)
            .sum()
    }
}

/// for (TR, R), the Mul operation is defined as each element of TR is
/// multiplied by R
impl<R: Ring, const K: usize> ops::Mul<R> for TR<R, K> {
    type Output = TR<R, K>;
    fn mul(self, other: R) -> TR<R, K> {
        Self(self.0.iter().map(|s| s.clone() * other.clone()).collect())
    }
}
impl<R: Ring, const K: usize> ops::Mul<&R> for &TR<R, K> {
    type Output = TR<R, K>;
    fn mul(self, other: &R) -> TR<R, K> {
        TR::<R, K>(self.0.iter().map(|s| s.clone() * other.clone()).collect())
    }
}
