use rand::{distributions::Distribution, Rng};
use std::fmt::Debug;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

/// Represents a ring element. Currently implemented by ring_nq.rs#Rq and
/// ring_torus.rs#Tn. Is not a 'pure algebraic ring', but more a custom trait
/// definition which includes methods like `mod_switch`.
// assumed to be mod (X^N +1)
pub trait Ring:
    Sized
    + Add<Output = Self>
    + AddAssign
    + Sum
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self> // internal product
    + Mul<u64, Output = Self> // scalar mul, external product
    + Mul<Self::C, Output = Self>
    + Neg<Output = Self>
    + PartialEq
    + Debug
    + Clone
    + Copy
    + Sum<<Self as Add>::Output>
    + Sum<<Self as Mul>::Output>
{
    /// C defines the coefficient type
    type C: Debug + Clone;

    const Q: u64;
    const N: usize;

    fn coeffs(&self) -> Vec<Self::C>;
    fn zero() -> Self;
    // note/wip/warning: dist (0,q) with f64, will output more '0=q' elements than other values
    fn rand(rng: impl Rng, dist: impl Distribution<f64>) -> Self;

    fn from_vec(coeffs: Vec<Self::C>) -> Self;

    fn decompose(&self, beta: u32, l: u32) -> Vec<Self>;

    fn remodule<const P: u64>(&self) -> impl Ring;
    fn mod_switch<const P: u64>(&self) -> impl Ring;

    /// returns [ [(num/den) * self].round() ] mod q
    /// ie. performs the multiplication and division over f64, and then it
    /// rounds the result, only applying the mod Q (if the ring is mod Q) at the
    /// end.
    fn mul_div_round(&self, num: u64, den: u64) -> Self;
}
