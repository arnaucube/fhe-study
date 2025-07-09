use anyhow::Result;
use rand::{distributions::Distribution, Rng};
use std::fmt::Debug;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

/// represents a ring element
pub trait Ring:
    Sized
    + Add<Output = Self>
    + AddAssign
    + Sum
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + Mul<u64, Output = Self> // scalar mul
    + PartialEq
    + Debug
    + Clone
    + Sum<<Self as Add>::Output>
    + Sum<<Self as Mul>::Output>
{
    /// C defines the coefficient type
    type C: Debug + Clone;

    fn coeffs(&self) -> Vec<Self::C>;
    fn zero() -> Self;
    fn rand(rng: impl Rng, dist: impl Distribution<f64>) -> Self;
    // note/wip/warning: dist (0,q) with f64, will output more '0=q' elements than other values
}
