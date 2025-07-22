//! Polynomial ring Z_q[X]/(X^N+1)
//!

use anyhow::{anyhow, Result};
use rand::{distributions::Distribution, Rng};
use std::array;
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

use crate::ntt::NTT;
use crate::zq::{modulus_u64, Zq};

use crate::Ring;

// NOTE: currently using fixed-size arrays, but pending to see if with
// real-world parameters the stack can keep up; if not will move everything to
// use Vec.
/// PolynomialRing element, where the PolynomialRing is R = Z_q[X]/(X^n +1)
/// The implementation assumes that q is prime.
#[derive(Clone, Copy)]
pub struct Rq<const Q: u64, const N: usize> {
    pub(crate) coeffs: [Zq<Q>; N],

    // evals are set when doig a PRxPR multiplication, so it can be reused in future
    // multiplications avoiding recomputing it
    pub(crate) evals: Option<[Zq<Q>; N]>,
}

impl<const Q: u64, const N: usize> Ring for Rq<Q, N> {
    type C = Zq<Q>;

    fn coeffs(&self) -> Vec<Self::C> {
        self.coeffs.to_vec()
    }
    fn zero() -> Self {
        let coeffs = array::from_fn(|_| Zq::zero());
        Self {
            coeffs,
            evals: None,
        }
    }
    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        // let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_u64(dist.sample(&mut rng)));
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Self::C::rand(&mut rng, &dist));
        Self {
            coeffs,
            evals: None,
        }
    }

    fn from_vec(coeffs: Vec<Zq<Q>>) -> Self {
        let mut p = coeffs;
        modulus::<Q, N>(&mut p);
        let coeffs = array::from_fn(|i| p[i]);
        Self {
            coeffs,
            evals: None,
        }
    }

    // returns the decomposition of each polynomial coefficient, such
    // decomposition will be a vecotor of length N, containint N vectors of Zq
    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        let elems: Vec<Vec<Zq<Q>>> = self.coeffs.iter().map(|r| r.decompose(beta, l)).collect();
        // transpose it
        let r: Vec<Vec<Zq<Q>>> = (0..elems[0].len())
            .map(|i| (0..elems.len()).map(|j| elems[j][i]).collect())
            .collect();
        // convert it to Rq<Q,N>
        r.iter().map(|a_i| Self::from_vec(a_i.clone())).collect()
    }

    // returns [ [(num/den) * self].round() ] mod q
    // ie. performs the multiplication and division over f64, and then it rounds the
    // result, only applying the mod Q at the end
    fn mul_div_round(&self, num: u64, den: u64) -> Self {
        let r: Vec<f64> = self
            .coeffs()
            .iter()
            .map(|e| ((num as f64 * e.0 as f64) / den as f64).round())
            .collect();
        Rq::<Q, N>::from_vec_f64(r)
    }
}

impl<const Q: u64, const N: usize> From<crate::ring_n::R<N>> for Rq<Q, N> {
    fn from(r: crate::ring_n::R<N>) -> Self {
        Self::from_vec(
            r.coeffs()
                .iter()
                .map(|e| Zq::<Q>::from_f64(*e as f64))
                .collect(),
        )
    }
}

// apply mod (X^N+1)
pub fn modulus<const Q: u64, const N: usize>(p: &mut Vec<Zq<Q>>) {
    if p.len() < N {
        return;
    }
    for i in N..p.len() {
        p[i - N] = p[i - N].clone() - p[i].clone();
        p[i] = Zq(0);
    }
    p.truncate(N);
}

// PR stands for PolynomialRing
impl<const Q: u64, const N: usize> Rq<Q, N> {
    pub fn coeffs(&self) -> [Zq<Q>; N] {
        self.coeffs
    }
    pub fn compute_evals(&mut self) {
        self.evals = Some(NTT::<Q, N>::ntt(self.coeffs));
    }
    pub fn to_r(self) -> crate::R<N> {
        crate::R::<N>::from(self)
    }

    // TODO rm since it is implemented in Ring trait impl
    // pub fn zero() -> Self {
    //     let coeffs = array::from_fn(|_| Zq::zero());
    //     Self {
    //         coeffs,
    //         evals: None,
    //     }
    // }
    // this method is mostly for tests
    pub fn from_vec_u64(coeffs: Vec<u64>) -> Self {
        let coeffs_mod_q = coeffs.iter().map(|c| Zq::from_u64(*c)).collect();
        Self::from_vec(coeffs_mod_q)
    }
    pub fn from_vec_f64(coeffs: Vec<f64>) -> Self {
        let coeffs_mod_q = coeffs.iter().map(|c| Zq::from_f64(*c)).collect();
        Self::from_vec(coeffs_mod_q)
    }
    pub fn from_vec_i64(coeffs: Vec<i64>) -> Self {
        let coeffs_mod_q = coeffs.iter().map(|c| Zq::from_f64(*c as f64)).collect();
        Self::from_vec(coeffs_mod_q)
    }
    pub fn new(coeffs: [Zq<Q>; N], evals: Option<[Zq<Q>; N]>) -> Self {
        Self { coeffs, evals }
    }

    pub fn rand_abs(mut rng: impl Rng, dist: impl Distribution<f64>) -> Result<Self> {
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_f64(dist.sample(&mut rng).abs()));
        Ok(Self {
            coeffs,
            evals: None,
        })
    }
    pub fn rand_f64_abs(mut rng: impl Rng, dist: impl Distribution<f64>) -> Result<Self> {
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_f64(dist.sample(&mut rng).abs()));
        Ok(Self {
            coeffs,
            evals: None,
        })
    }
    pub fn rand_f64(mut rng: impl Rng, dist: impl Distribution<f64>) -> Result<Self> {
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_f64(dist.sample(&mut rng)));
        Ok(Self {
            coeffs,
            evals: None,
        })
    }
    pub fn rand_u64(mut rng: impl Rng, dist: impl Distribution<u64>) -> Result<Self> {
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_u64(dist.sample(&mut rng)));
        Ok(Self {
            coeffs,
            evals: None,
        })
    }
    // WIP. returns random v \in {0,1}. // TODO {-1, 0, 1}
    pub fn rand_bin(mut rng: impl Rng, dist: impl Distribution<bool>) -> Result<Self> {
        let coeffs: [Zq<Q>; N] = array::from_fn(|_| Zq::from_bool(dist.sample(&mut rng)));
        Ok(Rq {
            coeffs,
            evals: None,
        })
    }
    // Warning: this method will behave differently depending on the values P and Q:
    // if Q<P, it just 'renames' the modulus parameter to P
    // if Q>=P, it crops to mod P
    pub fn remodule<const P: u64>(&self) -> Rq<P, N> {
        Rq::<P, N>::from_vec_u64(self.coeffs().iter().map(|m_i| m_i.0).collect())
    }

    /// perform the mod switch operation from Q to Q', where Q2=Q'
    pub fn mod_switch<const Q2: u64>(&self) -> Rq<Q2, N> {
        Rq::<Q2, N> {
            coeffs: array::from_fn(|i| self.coeffs[i].mod_switch::<Q2>()),
            evals: None,
        }
    }

    // applies mod(T) to all coefficients of self
    pub fn coeffs_mod<const T: u64>(&self) -> Self {
        Rq::<Q, N>::from_vec_u64(
            self.coeffs()
                .iter()
                .map(|m_i| modulus_u64::<T>(m_i.0))
                .collect(),
        )
    }

    // TODO review if needed, or if with this interface
    pub fn mul_by_matrix(&self, m: &Vec<Vec<Zq<Q>>>) -> Result<Vec<Zq<Q>>> {
        matrix_vec_product(m, &self.coeffs.to_vec())
    }
    pub fn mul_by_zq(&self, s: &Zq<Q>) -> Self {
        Self {
            coeffs: array::from_fn(|i| self.coeffs[i] * *s),
            evals: None,
        }
    }
    pub fn mul_by_u64(&self, s: u64) -> Self {
        let s = Zq::from_u64(s);
        Self {
            coeffs: array::from_fn(|i| self.coeffs[i] * s),
            // coeffs: self.coeffs.iter().map(|&e| e * s).collect(),
            evals: None,
        }
    }
    pub fn mul_by_f64(&self, s: f64) -> Self {
        Self {
            coeffs: array::from_fn(|i| Zq::from_f64(self.coeffs[i].0 as f64 * s)),
            evals: None,
        }
    }

    pub fn mul(&mut self, rhs: &mut Self) -> Self {
        mul_mut(self, rhs)
    }
    // divides by the given scalar 's' and rounds, returning a Rq<Q,N>
    // TODO rm
    pub fn div_round(&self, s: u64) -> Self {
        let r: Vec<f64> = self
            .coeffs()
            .iter()
            .map(|e| (e.0 as f64 / s as f64).round())
            .collect();
        Rq::<Q, N>::from_vec_f64(r)
    }

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // TODO simplify
        let mut str = "";
        let mut zero = true;
        for (i, coeff) in self.coeffs.iter().enumerate().rev() {
            if coeff.0 == 0 {
                continue;
            }
            zero = false;
            f.write_str(str)?;
            if coeff.0 != 1 {
                f.write_str(coeff.0.to_string().as_str())?;
                if i > 0 {
                    f.write_str("*")?;
                }
            }
            if coeff.0 == 1 && i == 0 {
                f.write_str(coeff.0.to_string().as_str())?;
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

        f.write_str(" mod Z_")?;
        f.write_str(Q.to_string().as_str())?;
        f.write_str("/(X^")?;
        f.write_str(N.to_string().as_str())?;
        f.write_str("+1)")?;
        Ok(())
    }

    pub fn infinity_norm(&self) -> u64 {
        self.coeffs()
            .iter()
            .map(|x| if x.0 > (Q / 2) { Q - x.0 } else { x.0 })
            .fold(0, |a, b| a.max(b))
    }
    pub fn mod_centered_q(&self) -> crate::ring_n::R<N> {
        self.to_r().mod_centered_q::<Q>()
    }
}
pub fn matrix_vec_product<const Q: u64>(m: &Vec<Vec<Zq<Q>>>, v: &Vec<Zq<Q>>) -> Result<Vec<Zq<Q>>> {
    // assert_eq!(m.len(), m[0].len()); // TODO change to returning err
    // assert_eq!(m.len(), v.len());
    if m.len() != m[0].len() {
        return Err(anyhow!("expected 'm' to be a square matrix"));
    }
    if m.len() != v.len() {
        return Err(anyhow!(
            "m.len: {} should be equal to v.len(): {}",
            m.len(),
            v.len(),
        ));
    }

    Ok(m.iter()
        .map(|row| {
            row.iter()
                .zip(v.iter())
                .map(|(&row_i, &v_i)| row_i * v_i)
                .sum()
        })
        .collect::<Vec<Zq<Q>>>())
}
pub fn transpose<const Q: u64>(m: &[Vec<Zq<Q>>]) -> Vec<Vec<Zq<Q>>> {
    // TODO case when m[0].len()=0
    // TODO non square matrix
    let mut r: Vec<Vec<Zq<Q>>> = vec![vec![Zq(0); m[0].len()]; m.len()];
    for (i, m_row) in m.iter().enumerate() {
        for (j, m_ij) in m_row.iter().enumerate() {
            r[j][i] = *m_ij;
        }
    }
    r
}

impl<const Q: u64, const N: usize> PartialEq for Rq<Q, N> {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}
impl<const Q: u64, const N: usize> Add<Rq<Q, N>> for Rq<Q, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            coeffs: array::from_fn(|i| self.coeffs[i] + rhs.coeffs[i]),
            evals: None,
        }
        // Self {
        //     coeffs: self
        //         .coeffs
        //         .iter()
        //         .zip(rhs.coeffs)
        //         .map(|(a, b)| *a + b)
        //         .collect(),
        //     evals: None,
        // }
        // Self(r.iter_mut().map(|e| e.r#mod()).collect()) // TODO mod should happen auto in +
    }
}
impl<const Q: u64, const N: usize> Add<&Rq<Q, N>> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn add(self, rhs: &Rq<Q, N>) -> Self::Output {
        Rq {
            coeffs: array::from_fn(|i| self.coeffs[i] + rhs.coeffs[i]),
            evals: None,
        }
    }
}
impl<const Q: u64, const N: usize> AddAssign for Rq<Q, N> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl<const Q: u64, const N: usize> Sum<Rq<Q, N>> for Rq<Q, N> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut acc = Rq::<Q, N>::zero();
        for e in iter {
            acc += e;
        }
        acc
    }
}

impl<const Q: u64, const N: usize> Sub<Rq<Q, N>> for Rq<Q, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            coeffs: array::from_fn(|i| self.coeffs[i] - rhs.coeffs[i]),
            evals: None,
        }
    }
}
impl<const Q: u64, const N: usize> Sub<&Rq<Q, N>> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn sub(self, rhs: &Rq<Q, N>) -> Self::Output {
        Rq {
            coeffs: array::from_fn(|i| self.coeffs[i] - rhs.coeffs[i]),
            evals: None,
        }
    }
}
impl<const Q: u64, const N: usize> SubAssign for Rq<Q, N> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl<const Q: u64, const N: usize> Mul<Rq<Q, N>> for Rq<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        mul(&self, &rhs)
    }
}
impl<const Q: u64, const N: usize> Mul<&Rq<Q, N>> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn mul(self, rhs: &Rq<Q, N>) -> Self::Output {
        mul(self, rhs)
    }
}

// mul by Zq element
impl<const Q: u64, const N: usize> Mul<Zq<Q>> for Rq<Q, N> {
    type Output = Self;

    fn mul(self, s: Zq<Q>) -> Self {
        self.mul_by_zq(&s)
    }
}
impl<const Q: u64, const N: usize> Mul<&Zq<Q>> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn mul(self, s: &Zq<Q>) -> Self::Output {
        self.mul_by_zq(s)
    }
}
// mul by u64
impl<const Q: u64, const N: usize> Mul<u64> for Rq<Q, N> {
    type Output = Self;

    fn mul(self, s: u64) -> Self {
        self.mul_by_u64(s)
    }
}
impl<const Q: u64, const N: usize> Mul<&u64> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn mul(self, s: &u64) -> Self::Output {
        self.mul_by_u64(*s)
    }
}
// mul by f64
impl<const Q: u64, const N: usize> Mul<f64> for Rq<Q, N> {
    type Output = Self;

    fn mul(self, s: f64) -> Self {
        self.mul_by_f64(s)
    }
}
impl<const Q: u64, const N: usize> Mul<&f64> for &Rq<Q, N> {
    type Output = Rq<Q, N>;

    fn mul(self, s: &f64) -> Self::Output {
        self.mul_by_f64(*s)
    }
}

impl<const Q: u64, const N: usize> Neg for Rq<Q, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            coeffs: array::from_fn(|i| -self.coeffs[i]),
            evals: None,
        }
    }
}

// note: this assumes that Q is prime
fn mul_mut<const Q: u64, const N: usize>(lhs: &mut Rq<Q, N>, rhs: &mut Rq<Q, N>) -> Rq<Q, N> {
    // reuse evaluations if already computed
    if !lhs.evals.is_some() {
        lhs.evals = Some(NTT::<Q, N>::ntt(lhs.coeffs));
    };
    if !rhs.evals.is_some() {
        rhs.evals = Some(NTT::<Q, N>::ntt(rhs.coeffs));
    };
    let lhs_evals = lhs.evals.unwrap();
    let rhs_evals = rhs.evals.unwrap();

    let c_ntt: [Zq<Q>; N] = array::from_fn(|i| lhs_evals[i] * rhs_evals[i]);
    let c = NTT::<Q, { N }>::intt(c_ntt);
    Rq::new(c, Some(c_ntt))
}
// note: this assumes that Q is prime
// TODO impl karatsuba for non-prime Q
fn mul<const Q: u64, const N: usize>(lhs: &Rq<Q, N>, rhs: &Rq<Q, N>) -> Rq<Q, N> {
    // reuse evaluations if already computed
    let lhs_evals = if lhs.evals.is_some() {
        lhs.evals.unwrap()
    } else {
        NTT::<Q, N>::ntt(lhs.coeffs)
    };
    let rhs_evals = if rhs.evals.is_some() {
        rhs.evals.unwrap()
    } else {
        NTT::<Q, N>::ntt(rhs.coeffs)
    };

    let c_ntt: [Zq<Q>; N] = array::from_fn(|i| lhs_evals[i] * rhs_evals[i]);
    let c = NTT::<Q, { N }>::intt(c_ntt);
    Rq::new(c, Some(c_ntt))
}

impl<const Q: u64, const N: usize> fmt::Display for Rq<Q, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}
impl<const Q: u64, const N: usize> fmt::Debug for Rq<Q, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_ring() {
        // the test values used are generated with SageMath
        const Q: u64 = 7;
        const N: usize = 3;

        // p = 1x + 2x^2 + 3x^3 + 4 x^4 + 5 x^5 in R=Z_q[X]/(X^n +1)
        let p = Rq::<Q, N>::from_vec_u64(vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        // try with coefficients bigger than Q
        let p = Rq::<Q, N>::from_vec_u64(vec![0u64, 1, Q + 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        // try with other ring
        let p = Rq::<7, 4>::from_vec_u64(vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "3*x^3 + 2*x^2 + 3*x + 3 mod Z_7/(X^4+1)");

        let p = Rq::<Q, N>::from_vec_u64(vec![0u64, 0, 0, 0, 4, 5]);
        assert_eq!(p.to_string(), "2*x^2 + 3*x mod Z_7/(X^3+1)");

        let p = Rq::<Q, N>::from_vec_u64(vec![5u64, 4, 5, 2, 1, 0]);
        assert_eq!(p.to_string(), "5*x^2 + 3*x + 3 mod Z_7/(X^3+1)");

        let a = Rq::<Q, N>::from_vec_u64(vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(a.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        let b = Rq::<Q, N>::from_vec_u64(vec![5u64, 4, 3, 2, 1, 0]);
        assert_eq!(b.to_string(), "3*x^2 + 3*x + 3 mod Z_7/(X^3+1)");

        // add
        assert_eq!((a.clone() + b.clone()).to_string(), "0 mod Z_7/(X^3+1)");
        assert_eq!((&a + &b).to_string(), "0 mod Z_7/(X^3+1)");
        // assert_eq!((a.0.clone() + b.0.clone()).to_string(), "[0, 0, 0]"); // TODO

        // sub
        assert_eq!(
            (a.clone() - b.clone()).to_string(),
            "x^2 + x + 1 mod Z_7/(X^3+1)"
        );
    }

    #[test]
    fn test_mul() -> Result<()> {
        const Q: u64 = 2u64.pow(16) + 1;
        const N: usize = 4;

        let a: [u64; N] = [1u64, 2, 3, 4];
        let b: [u64; N] = [1u64, 2, 3, 4];
        let c: [u64; N] = [65513, 65517, 65531, 20];
        test_mul_opt::<Q, N>(a, b, c)?;

        let a: [u64; N] = [0u64, 0, 0, 2];
        let b: [u64; N] = [0u64, 0, 0, 2];
        let c: [u64; N] = [0u64, 0, 65533, 0];
        test_mul_opt::<Q, N>(a, b, c)?;

        // TODO more testvectors

        Ok(())
    }
    fn test_mul_opt<const Q: u64, const N: usize>(
        a: [u64; N],
        b: [u64; N],
        expected_c: [u64; N],
    ) -> Result<()> {
        let a: [Zq<Q>; N] = array::from_fn(|i| Zq::from_u64(a[i]));
        let mut a = Rq::new(a, None);
        let b: [Zq<Q>; N] = array::from_fn(|i| Zq::from_u64(b[i]));
        let mut b = Rq::new(b, None);
        let expected_c: [Zq<Q>; N] = array::from_fn(|i| Zq::from_u64(expected_c[i]));
        let expected_c = Rq::new(expected_c, None);

        let c = mul_mut(&mut a, &mut b);
        assert_eq!(c, expected_c);
        Ok(())
    }

    #[test]
    fn test_rq_decompose() -> Result<()> {
        const Q: u64 = 16;
        const N: usize = 4;
        let beta = 4;
        let l = 2;

        let a = Rq::<Q, N>::from_vec_u64(vec![7u64, 14, 3, 6]);
        let d = a.decompose(beta, l);

        assert_eq!(
            d[0].coeffs().to_vec(),
            vec![1u64, 3, 0, 1]
                .iter()
                .map(|e| Zq::<Q>::from_u64(*e))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            d[1].coeffs().to_vec(),
            vec![3u64, 2, 3, 2]
                .iter()
                .map(|e| Zq::<Q>::from_u64(*e))
                .collect::<Vec<_>>()
        );
        Ok(())
    }
}
