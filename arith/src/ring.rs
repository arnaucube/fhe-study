//! Polynomial ring Z[X]/(X^N+1)
//!

use std::array;
use std::fmt;
use std::ops;

// PolynomialRing element, where the PolynomialRing is R = Z[X]/(X^n +1)
#[derive(Clone, Copy)]
pub struct R<const N: usize>(pub [i64; N]);

impl<const Q: u64, const N: usize> From<crate::ringq::Rq<Q, N>> for R<N> {
    fn from(rq: crate::ringq::Rq<Q, N>) -> Self {
        Self::from_vec_u64(rq.coeffs().to_vec().iter().map(|e| e.0).collect())
    }
}

impl<const N: usize> R<N> {
    pub fn coeffs(&self) -> [i64; N] {
        self.0
    }
    pub fn to_rq<const Q: u64>(self) -> crate::Rq<Q, N> {
        crate::Rq::<Q, N>::from(self)
    }

    pub fn from_vec(coeffs: Vec<i64>) -> Self {
        let mut p = coeffs;
        modulus::<N>(&mut p);
        Self(array::from_fn(|i| p[i]))
    }
    // this method is mostly for tests
    pub fn from_vec_u64(coeffs: Vec<u64>) -> Self {
        let coeffs_i64 = coeffs.iter().map(|c| *c as i64).collect();
        Self::from_vec(coeffs_i64)
    }
    pub fn from_vec_f64(coeffs: Vec<f64>) -> Self {
        let coeffs_i64 = coeffs.iter().map(|c| c.round() as i64).collect();
        Self::from_vec(coeffs_i64)
    }
    pub fn new(coeffs: [i64; N]) -> Self {
        Self(coeffs)
    }
    pub fn mul_by_i64(&self, s: i64) -> Self {
        Self(array::from_fn(|i| self.0[i] * s))
    }
    // performs the multiplication and division over f64, and then it rounds the
    // result, only applying the mod Q at the end
    pub fn mul_div_round<const Q: u64>(&self, num: u64, den: u64) -> crate::Rq<Q, N> {
        let r: Vec<f64> = self
            .coeffs()
            .iter()
            .map(|e| ((num as f64 * *e as f64) / den as f64).round())
            .collect();
        crate::Rq::<Q, N>::from_vec_f64(r)
    }
}

pub fn mul_div_round<const Q: u64, const N: usize>(
    v: Vec<i64>,
    num: u64,
    den: u64,
) -> crate::Rq<Q, N> {
    // dbg!(&v);
    let r: Vec<f64> = v
        .iter()
        .map(|e| ((num as f64 * *e as f64) / den as f64).round())
        .collect();
    // dbg!(&r);
    crate::Rq::<Q, N>::from_vec_f64(r)
}

// TODO rename to make it clear that is not mod q, but mod X^N+1
// apply mod (X^N+1)
pub fn modulus<const N: usize>(p: &mut Vec<i64>) {
    if p.len() < N {
        return;
    }
    for i in N..p.len() {
        p[i - N] = p[i - N].clone() - p[i].clone();
        p[i] = 0;
    }
    p.truncate(N);
}
pub fn modulus_i128<const N: usize>(p: &mut Vec<i128>) {
    if p.len() < N {
        return;
    }
    for i in N..p.len() {
        p[i - N] = p[i - N].clone() - p[i].clone();
        p[i] = 0;
    }
    p.truncate(N);
}

impl<const N: usize> PartialEq for R<N> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<const N: usize> ops::Add<R<N>> for R<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}
impl<const N: usize> ops::Add<&R<N>> for &R<N> {
    type Output = R<N>;

    fn add(self, rhs: &R<N>) -> Self::Output {
        R(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}
impl<const N: usize> ops::Sub<R<N>> for R<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}
impl<const N: usize> ops::Sub<&R<N>> for &R<N> {
    type Output = R<N>;

    fn sub(self, rhs: &R<N>) -> Self::Output {
        R(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}
impl<const N: usize> ops::Mul<R<N>> for R<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        naive_poly_mul(&self, &rhs)
    }
}
impl<const N: usize> ops::Mul<&R<N>> for &R<N> {
    type Output = R<N>;

    fn mul(self, rhs: &R<N>) -> Self::Output {
        naive_poly_mul(self, rhs)
    }
}

// TODO WIP
pub fn naive_poly_mul<const N: usize>(poly1: &R<N>, poly2: &R<N>) -> R<N> {
    let poly1: Vec<i128> = poly1.0.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.0.iter().map(|c| *c as i128).collect();
    let mut result: Vec<i128> = vec![0; (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    // R::<N>::from_vec(result.iter().map(|c| *c as i64).collect())
    modulus_i128::<N>(&mut result);
    // dbg!(&result);
    // dbg!(R::<N>(array::from_fn(|i| result[i] as i64)).coeffs());

    // sanity check: check that there are no coeffs > i64_max
    assert_eq!(
        result,
        R::<N>(array::from_fn(|i| result[i] as i64))
            .coeffs()
            .iter()
            .map(|c| *c as i128)
            .collect::<Vec<_>>()
    );
    R(array::from_fn(|i| result[i] as i64))
}
pub fn naive_mul_2<const N: usize>(poly1: &Vec<i128>, poly2: &Vec<i128>) -> Vec<i128> {
    let mut result: Vec<i128> = vec![0; (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // apply mod (X^N + 1))
    // R::<N>::from_vec(result.iter().map(|c| *c as i64).collect())
    modulus_i128::<N>(&mut result);
    result
}

pub fn naive_mul<const N: usize>(poly1: &R<N>, poly2: &R<N>) -> Vec<i64> {
    let poly1: Vec<i128> = poly1.0.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.0.iter().map(|c| *c as i128).collect();
    let mut result = vec![0; (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }
    result.iter().map(|c| *c as i64).collect()
}
pub fn naive_mul_TMP<const N: usize>(poly1: &R<N>, poly2: &R<N>) -> Vec<i64> {
    let poly1: Vec<i128> = poly1.0.iter().map(|c| *c as i128).collect();
    let poly2: Vec<i128> = poly2.0.iter().map(|c| *c as i128).collect();
    let mut result: Vec<i128> = vec![0; (N * 2) - 1];
    for i in 0..N {
        for j in 0..N {
            result[i + j] = result[i + j] + poly1[i] * poly2[j];
        }
    }

    // dbg!(&result);
    modulus_i128::<N>(&mut result);
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
pub fn mod_centered_q<const Q: u64, const N: usize>(p: Vec<i128>) -> R<N> {
    let q: i128 = Q as i128;
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
    R::<N>::from_vec(r.iter().map(|v| *v as i64).collect::<Vec<i64>>())
}

// mul by u64
impl<const N: usize> ops::Mul<u64> for R<N> {
    type Output = Self;

    fn mul(self, s: u64) -> Self {
        self.mul_by_i64(s as i64)
    }
}
impl<const N: usize> ops::Mul<&u64> for &R<N> {
    type Output = R<N>;

    fn mul(self, s: &u64) -> Self::Output {
        self.mul_by_i64(*s as i64)
    }
}

impl<const N: usize> ops::Neg for R<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(array::from_fn(|i| -self.0[i]))
    }
}

impl<const N: usize> R<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut str = "";
        let mut zero = true;
        for (i, coeff) in self.0.iter().enumerate().rev() {
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
        f.write_str(N.to_string().as_str())?;
        f.write_str("+1)")?;
        Ok(())
    }
}
impl<const N: usize> fmt::Display for R<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}
impl<const N: usize> fmt::Debug for R<N> {
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
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 2;
        let q: i64 = Q as i64;

        // test vectors generated with SageMath
        let a: [i64; N] = [q - 1, q - 1];
        let b: [i64; N] = [q - 1, q - 1];
        let c: [i64; N] = [0, 8589934592];
        test_mul_opt::<Q, N>(a, b, c)?;

        let a: [i64; N] = [1, q - 1];
        let b: [i64; N] = [1, q - 1];
        let c: [i64; N] = [-4294967295, 131072];
        test_mul_opt::<Q, N>(a, b, c)?;

        Ok(())
    }
    fn test_mul_opt<const Q: u64, const N: usize>(
        a: [i64; N],
        b: [i64; N],
        expected_c: [i64; N],
    ) -> Result<()> {
        let mut a = R::new(a);
        let mut b = R::new(b);
        dbg!(&a);
        dbg!(&b);
        let expected_c = R::new(expected_c);

        let mut c = naive_mul(&mut a, &mut b);
        modulus::<N>(&mut c);
        dbg!(R::<N>::from_vec(c.clone()));
        assert_eq!(c, expected_c.0.to_vec());
        Ok(())
    }
}
