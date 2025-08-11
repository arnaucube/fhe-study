//! Polynomial ring Z[X]/(X^N+1)
//!

use anyhow::Result;
use itertools::zip_eq;
use rand::{distributions::Distribution, Rng};
use std::array;
use std::borrow::Borrow;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::Ring;

// TODO rename to not have name conflicts with the Ring trait (R: Ring)
// PolynomialRing element, where the PolynomialRing is R = Z[X]/(X^n +1)
#[derive(Clone)]
pub struct R {
    pub n: usize,
    pub coeffs: Vec<i64>,
}

// impl<const N: usize> Ring for R<N> {
impl R {
    // type C = i64;
    // type Params = usize; // n

    // const Q: u64 = i64::MAX as u64; // WIP
    // const N: usize = N;

    pub fn coeffs(&self) -> Vec<i64> {
        self.coeffs.clone()
    }
    fn zero(n: usize) -> Self {
        Self {
            n,
            coeffs: vec![0i64; n],
        }
    }
    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>, n: usize) -> Self {
        // let coeffs: [i64; N] = array::from_fn(|_| Self::C::rand(&mut rng, &dist));
        // let coeffs: [i64; N] = array::from_fn(|_| dist.sample(&mut rng).round() as i64);
        Self {
            n,
            coeffs: std::iter::repeat_with(|| dist.sample(&mut rng).round() as i64)
                .take(n)
                .collect(),
        }
        // let coeffs: [C; N] = array::from_fn(|_| Zq::from_u64(dist.sample(&mut rng)));
        // Self(coeffs)
    }

    pub fn from_vec(n: usize, coeffs: Vec<i64>) -> Self {
        let mut p = coeffs;
        modulus(n, &mut p);
        Self { n, coeffs: p }
    }

    /*
    // returns the decomposition of each polynomial coefficient
    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        unimplemented!();
        // array::from_fn(|i| self.coeffs[i].decompose(beta, l))
    }

    fn remodule<const P: u64>(&self) -> impl Ring {
        unimplemented!()
    }
    fn mod_switch<const P: u64, const M: usize>(&self) -> R<N> {
        unimplemented!()
    }

    /// performs the multiplication and division over f64, and then it rounds the
    /// result, only applying the mod Q at the end
    fn mul_div_round(&self, num: u64, den: u64) -> Self {
        unimplemented!()
        // fn mul_div_round<const Q: u64>(&self, num: u64, den: u64) -> crate::Rq<Q, N> {
        // let r: Vec<f64> = self
        //     .coeffs()
        //     .iter()
        //     .map(|e| ((num as f64 * *e as f64) / den as f64).round())
        //     .collect();
        // crate::Rq::<Q, N>::from_vec_f64(r)
    }
    */
}

impl From<crate::ring_nq::Rq> for R {
    fn from(rq: crate::ring_nq::Rq) -> Self {
        Self::from_vec_u64(rq.n, rq.coeffs().to_vec().iter().map(|e| e.v).collect())
    }
}

impl R {
    // pub fn coeffs(&self) -> [i64; N] {
    //     self.0
    // }
    pub fn to_rq(self, q: u64) -> crate::Rq {
        crate::Rq::from((q, self))
    }

    // this method is mostly for tests
    pub fn from_vec_u64(n: usize, coeffs: Vec<u64>) -> Self {
        let coeffs_i64 = coeffs.iter().map(|c| *c as i64).collect();
        Self::from_vec(n, coeffs_i64)
    }
    pub fn from_vec_f64(n: usize, coeffs: Vec<f64>) -> Self {
        let coeffs_i64 = coeffs.iter().map(|c| c.round() as i64).collect();
        Self::from_vec(n, coeffs_i64)
    }
    pub fn new(n: usize, coeffs: Vec<i64>) -> Self {
        assert_eq!(n, coeffs.len());
        Self { n, coeffs }
    }
    pub fn mul_by_i64(&self, s: i64) -> Self {
        Self {
            n: self.n,
            coeffs: self.coeffs.iter().map(|c_i| c_i * s).collect(),
        }
    }

    pub fn infinity_norm(&self) -> u64 {
        self.coeffs()
            .iter()
            // .map(|x| if x.0 > (Q / 2) { Q - x.0 } else { x.0 })
            .map(|x| x.abs() as u64)
            .fold(0, |a, b| a.max(b))
    }
    pub fn mod_centered_q(&self, q: u64) -> R {
        let q = q as i64;
        let r = self
            .coeffs
            .iter()
            .map(|v| {
                let mut res = v % q;
                if res > q / 2 {
                    res = res - q;
                }
                res
            })
            .collect::<Vec<i64>>();
        R::from_vec(self.n, r)
    }
}

pub fn mul_div_round(q: u64, n: usize, v: Vec<i64>, num: u64, den: u64) -> crate::Rq {
    // dbg!(&v);
    let r: Vec<f64> = v
        .iter()
        .map(|e| ((num as f64 * *e as f64) / den as f64).round())
        .collect();
    // dbg!(&r);
    crate::Rq::from_vec_f64(q, n, r)
}

// TODO rename to make it clear that is not mod q, but mod X^N+1
// apply mod (X^N+1)
pub fn modulus(n: usize, p: &mut Vec<i64>) {
    if p.len() < n {
        return;
    }
    for i in n..p.len() {
        p[i - n] = p[i - n].clone() - p[i].clone();
        p[i] = 0;
    }
    p.truncate(n);
}
pub fn modulus_i128(n: usize, p: &mut Vec<i128>) {
    if p.len() < n {
        return;
    }
    for i in n..p.len() {
        p[i - n] = p[i - n].clone() - p[i].clone();
        p[i] = 0;
    }
    p.truncate(n);
}

impl PartialEq for R {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs && self.n == other.n
    }
}
impl Add<R> for R {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.n, rhs.n);
        Self {
            n: self.n,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l + r)
                .collect(),
        }
    }
}
impl Add<&R> for &R {
    type Output = R;

    fn add(self, rhs: &R) -> Self::Output {
        // R(array::from_fn(|i| self.0[i] + rhs.0[i]))
        assert_eq!(self.n, rhs.n);
        R {
            n: self.n,
            coeffs: zip_eq(self.coeffs.clone(), rhs.coeffs.clone())
                .map(|(l, r)| l + r)
                .collect(),
        }
    }
}
impl AddAssign for R {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.n, rhs.n);
        for i in 0..self.n {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl Sum<R> for R {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        // let mut acc = R::zero();
        // for e in iter {
        //     acc += e;
        // }
        // acc
        let first = iter.next().unwrap();
        iter.fold(first, |acc, x| acc + x)
    }
}

impl Sub<R> for R {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        // Self(array::from_fn(|i| self.0[i] - rhs.0[i]))
        assert_eq!(self.n, rhs.n);
        Self {
            n: self.n,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l - r)
                .collect(),
        }
    }
}
impl Sub<&R> for &R {
    type Output = R;

    fn sub(self, rhs: &R) -> Self::Output {
        // R(array::from_fn(|i| self.0[i] - rhs.0[i]))
        assert_eq!(self.n, rhs.n);
        R {
            n: self.n,
            coeffs: zip_eq(&self.coeffs, &rhs.coeffs)
                .map(|(l, r)| l - r)
                .collect(),
        }
    }
}

impl SubAssign for R {
    fn sub_assign(&mut self, rhs: Self) {
        assert_eq!(self.n, rhs.n);
        for i in 0..self.n {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl Mul<R> for R {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        naive_poly_mul(&self, &rhs)
    }
}
impl Mul<&R> for &R {
    type Output = R;

    fn mul(self, rhs: &R) -> Self::Output {
        naive_poly_mul(self, rhs)
    }
}

// TODO WIP
pub fn naive_poly_mul(poly1: &R, poly2: &R) -> R {
    assert_eq!(poly1.n, poly2.n);
    let n = poly1.n;

    let poly1: Vec<i128> = poly1.coeffs.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.coeffs.iter().map(|c| *c as i128).collect();
    let mut result: Vec<i128> = vec![0; (n * 2) - 1];
    for i in 0..n {
        for j in 0..n {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    // R::<N>::from_vec(result.iter().map(|c| *c as i64).collect())
    modulus_i128(n, &mut result);
    // dbg!(&result);
    // dbg!(R::<N>(array::from_fn(|i| result[i] as i64)).coeffs());

    let result_i64: Vec<i64> = result.iter().map(|c_i| *c_i as i64).collect();
    let r = R::from_vec(n, result_i64);
    // sanity check: check that there are no coeffs > i64_max
    assert_eq!(
        result,
        r.coeffs.iter().map(|c| *c as i128).collect::<Vec<_>>()
    );
    r
}
pub fn naive_mul_2(n: usize, poly1: &Vec<i128>, poly2: &Vec<i128>) -> Vec<i128> {
    let mut result: Vec<i128> = vec![0; (n * 2) - 1];
    for i in 0..n {
        for j in 0..n {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    // R::<N>::from_vec(result.iter().map(|c| *c as i64).collect())
    modulus_i128(n, &mut result);
    result
}

pub fn naive_mul(poly1: &R, poly2: &R) -> Vec<i64> {
    assert_eq!(poly1.n, poly2.n);
    let n = poly1.n;

    let poly1: Vec<i128> = poly1.coeffs.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.coeffs.iter().map(|c| *c as i128).collect();
    let mut result = vec![0; (n * 2) - 1];
    for i in 0..n {
        for j in 0..n {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }
    result.iter().map(|c| *c as i64).collect()
}
pub fn naive_mul_TMP(poly1: &R, poly2: &R) -> Vec<i64> {
    assert_eq!(poly1.n, poly2.n);
    let n = poly1.n;

    let poly1: Vec<i128> = poly1.coeffs.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.coeffs.iter().map(|c| *c as i128).collect();
    let mut result: Vec<i128> = vec![0; (n * 2) - 1];
    for i in 0..n {
        for j in 0..n {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // dbg!(&result);
    modulus_i128(n, &mut result);
    // for c_i in result.iter() {
    //     println!("---");
    //     println!("{:?}", &c_i);
    //     println!("{:?}", *c_i as i64);
    //     println!("{:?}", (*c_i as i64) as i128);
    //     assert_eq!(*c_i, (*c_i as i64) as i128, "{:?}", c_i);
    // }
    result.iter().map(|c| *c as i64).collect()
}

// wip
pub fn mod_centered_q(q: u64, n: usize, p: Vec<i128>) -> R {
    let q: i128 = q as i128;
    let r = p
        .iter()
        .map(|v| {
            let mut res = v % q;
            if res > q / 2 {
                res = res - q;
            }
            res
        })
        .collect::<Vec<i128>>();
    R::from_vec(n, r.iter().map(|v| *v as i64).collect::<Vec<i64>>())
}

impl Mul<i64> for R {
    type Output = Self;

    fn mul(self, s: i64) -> Self {
        self.mul_by_i64(s)
    }
}
// mul by u64
impl Mul<u64> for R {
    type Output = Self;

    fn mul(self, s: u64) -> Self {
        self.mul_by_i64(s as i64)
    }
}
impl Mul<&u64> for &R {
    type Output = R;

    fn mul(self, s: &u64) -> Self::Output {
        self.mul_by_i64(*s as i64)
    }
}

impl Neg for R {
    type Output = Self;

    fn neg(self) -> Self::Output {
        // Self(array::from_fn(|i| -self.0[i]))
        Self {
            n: self.n,
            coeffs: self.coeffs.iter().map(|c_i| -c_i).collect(),
        }
    }
}

impl R {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut str = "";
        let mut zero = true;
        for (i, coeff) in self.coeffs.iter().enumerate().rev() {
            if *coeff == 0 {
                continue;
            }
            zero = false;
            f.write_str(str)?;
            if *coeff != 1 {
                f.write_str(coeff.to_string().as_str())?;
                if i > 0 {
                    f.write_str("*")?;
                }
            }
            if *coeff == 1 && i == 0 {
                f.write_str(coeff.to_string().as_str())?;
            }
            if i == 1 {
                f.write_str("x")?;
            } else if i > 1 {
                f.write_str("x^")?;
                f.write_str(i.to_string().as_str())?;
            }
            str = " + ";
        }
        if zero {
            f.write_str("0")?;
        }

        f.write_str(" mod Z")?;
        f.write_str("/(X^")?;
        f.write_str(self.n.to_string().as_str())?;
        f.write_str("+1)")?;
        Ok(())
    }
}
impl fmt::Display for R {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}
impl fmt::Debug for R {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_mul() -> Result<()> {
        let n: usize = 2;
        let q: i64 = (2u64.pow(16) + 1) as i64;

        // test vectors generated with SageMath
        let a: Vec<i64> = vec![q - 1, q - 1];
        let b: Vec<i64> = vec![q - 1, q - 1];
        let c: Vec<i64> = vec![0, 8589934592];
        test_mul_opt(n, a, b, c)?;

        let a: Vec<i64> = vec![1, q - 1];
        let b: Vec<i64> = vec![1, q - 1];
        let c: Vec<i64> = vec![-4294967295, 131072];
        test_mul_opt(n, a, b, c)?;

        Ok(())
    }
    fn test_mul_opt(n: usize, a: Vec<i64>, b: Vec<i64>, expected_c: Vec<i64>) -> Result<()> {
        let mut a = R::new(n, a);
        let mut b = R::new(n, b);
        dbg!(&a);
        dbg!(&b);
        let expected_c = R::new(n, expected_c);

        let mut c = naive_mul(&mut a, &mut b);
        modulus(n, &mut c);
        dbg!(R::from_vec(n, c.clone()));
        assert_eq!(c, expected_c.coeffs);
        Ok(())
    }
}
