//! Complex

use rand::Rng;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct C<T> {
    pub re: T,
    pub im: T,
}

impl From<f64> for C<f64> {
    fn from(v: f64) -> Self {
        Self { re: v, im: 0_f64 }
    }
}

impl<T> C<T>
where
    T: Default + From<i32>,
{
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    pub fn rand(mut rng: impl Rng, max: u64) -> Self {
        Self::new(
            T::from(rng.gen_range(0..max) as i32),
            T::from(rng.gen_range(0..max) as i32),
        )
    }

    pub fn zero() -> C<T> {
        Self {
            re: T::from(0),
            im: T::from(0),
        }
    }
    pub fn one() -> C<T> {
        Self {
            re: T::from(1),
            im: T::from(0),
        }
    }
    pub fn i() -> C<T> {
        Self {
            re: T::from(0),
            im: T::from(1),
        }
    }
}

impl C<f64> {
    // cos & sin from Taylor series approximation, details at
    // https://en.wikipedia.org/wiki/Sine_and_cosine#Series_and_polynomials
    fn cos(x: f64) -> f64 {
        let mut r = 1.0;
        let mut term = 1.0;
        let mut n = 1;
        for _ in 0..10 {
            term *= -(x * x) / ((2 * n - 1) * (2 * n)) as f64;
            r += term;
            n += 1;
        }

        r
    }
    fn sin(x: f64) -> f64 {
        let mut r = x;
        let mut term = x;
        let mut n = 1;

        for _ in 0..10 {
            term *= -(x * x) / ((2 * n) * (2 * n + 1)) as f64;
            r += term;
            n += 1;
        }

        r
    }

    // e^(self))
    pub fn exp(self) -> Self {
        Self {
            re: Self::cos(self.im), // TODO WIP review
            im: Self::sin(self.im),
        }
    }
    pub fn pow(self, k: u32) -> Self {
        let mut k = k;
        if k == 0 {
            return Self::one();
        }
        let mut base = self.clone();
        while k & 1 == 0 {
            base = base.clone() * base;
            k >>= 1;
        }

        if k == 1 {
            return base;
        }

        let mut acc = base.clone();
        while k > 1 {
            k >>= 1;
            base = base.clone() * base;
            if k & 1 == 1 {
                acc = acc * base.clone();
            }
        }
        acc
    }
    pub fn modulus<const Q: u64>(self) -> Self {
        let q: f64 = Q as f64;
        let re = (self.re % q + q) % q;
        let im = (self.im % q + q) % q;
        Self { re, im }
    }
    pub fn modulus_centered<const Q: u64>(self) -> Self {
        let re = modulus_centered_f64::<Q>(self.re);
        let im = modulus_centered_f64::<Q>(self.im);
        Self { re, im }
    }
}

fn modulus_centered_f64<const Q: u64>(v: f64) -> f64 {
    let q = Q as f64;
    let mut res = v % q;
    if res > q / 2.0 {
        res = res - q;
    }
    res
}

impl<T: Default> Default for C<T> {
    fn default() -> Self {
        C {
            re: T::default(),
            im: T::default(),
        }
    }
}

impl<T> Add for C<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        C {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl<T> Sub for C<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        C {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl<T> Mul for C<T>
where
    T: Mul<Output = T> + Sub<Output = T> + Add<Output = T> + Copy,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        C {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl<T> Neg for C<T>
where
    T: Neg<Output = T> + Copy,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        C {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl<T> Div for C<T>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Neg<Output = T>,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        // (a+ib)/(c+id) = (ac + bd)/(c^2 + d^2) + i* (bc -ad)/(c^2 + d^2)
        let den = rhs.re * rhs.re + rhs.im * rhs.im;
        C {
            re: (self.re * rhs.re + self.im * rhs.im) / den,
            im: (self.im * rhs.re - self.re * rhs.im) / den,
        }
    }
}
impl<T> C<T>
where
    T: Neg<Output = T> + Copy,
{
    pub fn conj(&self) -> Self {
        C {
            re: self.re,
            im: -self.im,
        }
    }
}

impl<T: Add<Output = T> + Default + Copy> std::iter::Sum for C<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(C::default(), |acc, x| acc + x)
    }
}

impl C<f64> {
    pub fn abs(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

// poly mul with complex coefficients
pub fn naive_poly_mul<const N: usize>(poly1: &Vec<C<f64>>, poly2: &Vec<C<f64>>) -> Vec<C<f64>> {
    let mut result: Vec<C<f64>> = vec![C::<f64>::zero(); (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    // R::<N>::from_vec(result.iter().map(|c| *c as i64).collect())
    // modulus_i128::<N>(&mut result);
    // dbg!(&result);
    // dbg!(R::<N>(array::from_fn(|i| result[i] as i64)).coeffs());

    // R(array::from_fn(|i| result[i] as i64))
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_i() {
        assert_eq!(C::i(), C::new(0.0, 1.0));
    }

    #[test]
    fn test_add() {
        let a = C::new(2.0, 3.0);
        let b = C::new(1.0, 4.0);
        assert_eq!(a + b, C::new(3.0, 7.0));
    }

    #[test]
    fn test_sub() {
        let a = C::new(5.0, 7.0);
        let b = C::new(2.0, 3.0);
        assert_eq!(a - b, C::new(3.0, 4.0));
    }

    #[test]
    fn test_mult() {
        let a = C::new(1.0, 2.0);
        let b = C::new(3.0, 4.0);
        assert_eq!(a * b, C::new(-5.0, 10.0));
    }

    #[test]
    fn test_div() {
        let a: C<f64> = C::new(1.0, 2.0);
        let b: C<f64> = C::new(3.0, 4.0);
        let r = a / b;
        let expected = C::new(0.44, 0.08);
        let epsilon = 1e-2;
        assert!((r.re - expected.re).abs() < epsilon);
        assert!((r.im - expected.im).abs() < epsilon);
    }

    #[test]
    fn test_conj() {
        let a = C::new(3.0, -4.0);
        assert_eq!(a.conj(), C::new(3.0, 4.0));
        assert_eq!(a.conj().conj(), a);
    }

    #[test]
    fn test_neg() {
        let a = C::new(1.0, -2.0);
        assert_eq!(-a, C::new(-1.0, 2.0));
    }

    #[test]
    fn test_abs() {
        let a = C::new(3.0, 4.0);
        assert!((a.abs() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_exp() {
        // let a = C::new(3.0, 4.0);
        let pi = C::<f64>::from(std::f64::consts::PI);
        let n = 4;
        let a = ((C::<f64>::from(2f64) * pi * C::<f64>::i()) / C::<f64>::new(n as f64, 0f64)).exp();
        dbg!(&a);
        assert_eq!(a.exp(), a.exp());
    }
}
