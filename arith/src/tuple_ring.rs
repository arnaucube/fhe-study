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

use crate::{Ring, RingParam};

/// Tuple of K Ring (Rq) elements. We use Vec<R> to allocate it in the heap,
/// since if using a fixed-size array it would overflow the stack.
#[derive(Clone, Debug)]
pub struct TR<R: Ring> {
    pub k: usize,
    pub r: Vec<R>,
}
// TODO rm pub from Vec<R>, so that TR can not be created from a Vec with
// invalid length, since it has to be created using the `new` method.

impl<R: Ring> TR<R> {
    pub fn new(k: usize, r: Vec<R>) -> Self {
        assert_eq!(r.len(), k);
        Self { k, r }
    }
    pub fn zero(k: usize, r_params: &RingParam) -> Self {
        Self {
            k,
            r: (0..k).into_iter().map(|_| R::zero(r_params)).collect(),
        }
    }
    pub fn rand(
        mut rng: impl Rng,
        dist: impl Distribution<f64>,
        k: usize,
        r_params: &RingParam,
    ) -> Self {
        Self {
            k,
            r: (0..k)
                .into_iter()
                .map(|_| R::rand(&mut rng, &dist, r_params))
                .collect(),
        }
    }
    // returns the decomposition of each polynomial element
    pub fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        unimplemented!()
        // self.0.iter().map(|r| r.decompose(beta, l)).collect() // this is Vec<Vec<Vec<R::C>>>
    }
}

impl TR<crate::torus::T64> {
    pub fn mod_switch(&self, q2: u64) -> TR<crate::torus::T64> {
        TR::<crate::torus::T64> {
            k: self.k,
            r: self.r.iter().map(|c_i| c_i.mod_switch(q2)).collect(),
        }
    }
    // pub fn mod_switch(&self, Q2: u64) -> TR<crate::torus::T64, K> {
    //     TR(self.0.iter().map(|c_i| c_i.mod_switch(Q2)).collect())
    // }
}
impl TR<crate::ring_torus::Tn> {
    pub fn left_rotate(&self, h: usize) -> Self {
        TR {
            k: self.k,
            r: self.r.iter().map(|c_i| c_i.left_rotate(h)).collect(),
        }
    }
}

impl<R: Ring> TR<R> {
    pub fn iter(&self) -> std::slice::Iter<R> {
        self.r.iter()
    }
}

impl<R: Ring> Add<TR<R>> for TR<R> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        debug_assert_eq!(self.k, other.k);

        Self {
            k: self.k,
            r: zip_eq(self.r, other.r)
                .map(|(s, o)| s + o)
                .collect::<Vec<_>>(),
        }
    }
}

impl<R: Ring> Sub<TR<R>> for TR<R> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        debug_assert_eq!(self.k, other.k);

        Self {
            k: self.k,
            r: zip_eq(self.r, other.r).map(|(s, o)| s - o).collect(),
        }
    }
}

impl<R: Ring> Neg for TR<R> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            k: self.k,
            r: self.r.iter().map(|e_i| -*e_i).collect(),
        }
    }
}

/// for (TR,TR), the Mul operation is defined as the dot product:
/// for A, B \in R^k, result = Σ A_i * B_i \in R
impl<R: Ring> Mul<TR<R>> for TR<R> {
    type Output = R;
    fn mul(self, other: Self) -> R {
        debug_assert_eq!(self.k, other.k);

        zip_eq(self.r, other.r).map(|(s, o)| s * o).sum()
    }
}
impl<R: Ring> Mul<&TR<R>> for &TR<R> {
    type Output = R;
    fn mul(self, other: &TR<R>) -> R {
        debug_assert_eq!(self.k, other.k);

        zip_eq(self.r.clone(), other.r.clone())
            .map(|(s, o)| s * o)
            .sum()
    }
}

/// for (TR, R), the Mul operation is defined as each element of TR is
/// multiplied by R
impl<R: Ring> Mul<R> for TR<R> {
    type Output = TR<R>;
    fn mul(self, other: R) -> TR<R> {
        Self {
            k: self.k,
            r: self.r.iter().map(|s| s.clone() * other.clone()).collect(),
        }
    }
}
impl<R: Ring> Mul<&R> for &TR<R> {
    type Output = TR<R>;
    fn mul(self, other: &R) -> TR<R> {
        TR::<R> {
            k: self.k,
            r: self.r.iter().map(|s| s.clone() * other.clone()).collect(),
        }
    }
}
