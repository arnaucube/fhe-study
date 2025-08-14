//! Polynomial ring Z_q[X]/(X^N+1)
//!

use anyhow::{anyhow, Result};
use itertools::zip_eq;
use rand::{distributions::Distribution, Rng};
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

use crate::ntt::NTT;
use crate::zq::{modulus_u64, Zq};

use crate::{Ring, RingParam};

/// PolynomialRing element, where the PolynomialRing is R = Z_q[X]/(X^n +1)
/// The implementation assumes that q is prime.
#[derive(Clone)]
pub struct Rq {
    pub param: RingParam,

    pub(crate) coeffs: Vec<Zq>,

    // evals are set when doig a PRxPR multiplication, so it can be reused in future
    // multiplications avoiding recomputing it
    pub(crate) evals: Option<Vec<Zq>>,
}

impl Ring for Rq {
    type C = Zq;

    fn param(&self) -> RingParam {
        self.param
    }
    fn coeffs(&self) -> Vec<Self::C> {
        self.coeffs.to_vec()
    }
    fn zero(param: &RingParam) -> Self {
        Self {
            param: param.clone(),
            coeffs: vec![Zq::zero(param.q); param.n],
            evals: None,
        }
    }
    fn rand(mut rng: impl Rng, dist: impl Distribution<f64>, param: &RingParam) -> Self {
        Self {
            param: param.clone(),
            coeffs: std::iter::repeat_with(|| Self::C::rand(&mut rng, &dist, param.q))
                .take(param.n)
                .collect(),
            evals: None,
        }
    }

    fn from_vec(param: &RingParam, coeffs: Vec<Zq>) -> Self {
        let mut p = coeffs;
        modulus(param.q, param.n, &mut p);
        Self {
            param: param.clone(),
            coeffs: p,
            evals: None,
        }
    }

    // returns the decomposition of each polynomial coefficient, such
    // decomposition will be a vector of length N, containing N vectors of Zq
    fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        let elems: Vec<Vec<Zq>> = self.coeffs.iter().map(|r| r.decompose(beta, l)).collect();
        // transpose it
        let r: Vec<Vec<Zq>> = (0..elems[0].len())
            .map(|i| (0..elems.len()).map(|j| elems[j][i]).collect())
            .collect();
        // convert it to Rq<Q,N>
        r.iter()
            .map(|a_i| Self::from_vec(&self.param, a_i.clone()))
            .collect()
    }

    // Warning: this method will behave differently depending on the values P and Q:
    // if Q<P, it just 'renames' the modulus parameter to P
    // if Q>=P, it crops to mod P
    fn remodule(&self, p: u64) -> Rq {
        let param = RingParam {
            q: p,
            n: self.param.n,
        };
        Rq::from_vec_u64(&param, self.coeffs().iter().map(|m_i| m_i.v).collect())
    }

    /// perform the mod switch operation from Q to Q', where Q2=Q'
    fn mod_switch(&self, p: u64) -> Rq {
        let param = RingParam {
            q: p,
            n: self.param.n,
        };
        Rq {
            param,
            coeffs: self.coeffs.iter().map(|c_i| c_i.mod_switch(p)).collect(),
            evals: None,
        }
    }

    /// returns [ [(num/den) * self].round() ] mod q
    /// ie. performs the multiplication and division over f64, and then it rounds the
    /// result, only applying the mod Q at the end
    fn mul_div_round(&self, num: u64, den: u64) -> Self {
        let r: Vec<f64> = self
            .coeffs()
            .iter()
            .map(|e| ((num as f64 * e.v as f64) / den as f64).round())
            .collect();
        Rq::from_vec_f64(&self.param, r)
    }
}

impl From<(u64, crate::ring_n::R)> for Rq {
    fn from(qr: (u64, crate::ring_n::R)) -> Self {
        let (q, r) = qr;
        assert_eq!(r.n, r.coeffs.len());

        Self::from_vec(
            &RingParam { q, n: r.n },
            r.coeffs()
                .iter()
                .map(|e| Zq::from_f64(q, *e as f64))
                .collect(),
        )
    }
}

// apply mod (X^N+1)
pub fn modulus(q: u64, n: usize, p: &mut Vec<Zq>) {
    if p.len() < n {
        return;
    }
    for i in n..p.len() {
        p[i - n] = p[i - n].clone() - p[i].clone();
        p[i] = Zq::zero(q);
    }
    p.truncate(n);
}

impl Rq {
    pub fn coeffs(&self) -> Vec<Zq> {
        self.coeffs.clone()
    }
    pub fn compute_evals(&mut self) {
        self.evals = Some(NTT::ntt(self).coeffs);
        // TODO improve, ntt returns Rq but here just needs Vec<Zq>
    }
    pub fn to_r(self) -> crate::R {
        crate::R::from(self)
    }

    // this method is mostly for tests
    pub fn from_vec_u64(param: &RingParam, coeffs: Vec<u64>) -> Self {
        let coeffs_mod_q: Vec<Zq> = coeffs.iter().map(|c| Zq::from_u64(param.q, *c)).collect();
        Self::from_vec(param, coeffs_mod_q)
    }
    pub fn from_vec_f64(param: &RingParam, coeffs: Vec<f64>) -> Self {
        let coeffs_mod_q: Vec<Zq> = coeffs.iter().map(|c| Zq::from_f64(param.q, *c)).collect();
        Self::from_vec(param, coeffs_mod_q)
    }
    pub fn from_vec_i64(param: &RingParam, coeffs: Vec<i64>) -> Self {
        let coeffs_mod_q: Vec<Zq> = coeffs
            .iter()
            .map(|c| Zq::from_f64(param.q, *c as f64))
            .collect();
        Self::from_vec(param, coeffs_mod_q)
    }
    pub fn new(param: &RingParam, coeffs: Vec<Zq>, evals: Option<Vec<Zq>>) -> Self {
        Self {
            param: *param,
            coeffs,
            evals,
        }
    }

    pub fn rand_abs(
        mut rng: impl Rng,
        dist: impl Distribution<f64>,
        param: &RingParam,
    ) -> Result<Self> {
        Ok(Self {
            param: *param,
            coeffs: std::iter::repeat_with(|| Zq::from_f64(param.q, dist.sample(&mut rng).abs()))
                .take(param.n)
                .collect(),
            evals: None,
        })
    }
    pub fn rand_f64_abs(
        mut rng: impl Rng,
        dist: impl Distribution<f64>,
        param: &RingParam,
    ) -> Result<Self> {
        Ok(Self {
            param: *param,
            coeffs: std::iter::repeat_with(|| Zq::from_f64(param.q, dist.sample(&mut rng).abs()))
                .take(param.n)
                .collect(),
            evals: None,
        })
    }
    pub fn rand_f64(
        mut rng: impl Rng,
        dist: impl Distribution<f64>,
        param: &RingParam,
    ) -> Result<Self> {
        Ok(Self {
            param: *param,
            coeffs: std::iter::repeat_with(|| Zq::from_f64(param.q, dist.sample(&mut rng)))
                .take(param.n)
                .collect(),
            evals: None,
        })
    }
    pub fn rand_u64(
        mut rng: impl Rng,
        dist: impl Distribution<u64>,
        param: &RingParam,
    ) -> Result<Self> {
        Ok(Self {
            param: *param,
            coeffs: std::iter::repeat_with(|| Zq::from_u64(param.q, dist.sample(&mut rng)))
                .take(param.n)
                .collect(),
            evals: None,
        })
    }
    // WIP. returns random v \in {0,1}. // TODO {-1, 0, 1}
    pub fn rand_bin(
        mut rng: impl Rng,
        dist: impl Distribution<bool>,
        param: &RingParam,
    ) -> Result<Self> {
        Ok(Rq {
            param: *param,
            coeffs: std::iter::repeat_with(|| Zq::from_bool(param.q, dist.sample(&mut rng)))
                .take(param.n)
                .collect(),
            evals: None,
        })
    }
    // Warning: this method will behave differently depending on the values P and Q:
    // if Q<P, it just 'renames' the modulus parameter to P
    // if Q>=P, it crops to mod P
    // pub fn remodule<const P: u64>(&self) -> Rq<P, N> {
    //     Rq::<P, N>::from_vec_u64(self.coeffs().iter().map(|m_i| m_i.0).collect())
    // }

    // applies mod(T) to all coefficients of self
    pub fn coeffs_mod(&self, param: &RingParam, t: u64) -> Self {
        Rq::from_vec_u64(
            param,
            self.coeffs()
                .iter()
                .map(|m_i| modulus_u64(t, m_i.v))
                .collect(),
        )
    }

    // TODO review if needed, or if with this interface
    pub fn mul_by_matrix(&self, m: &Vec<Vec<Zq>>) -> Result<Vec<Zq>> {
        matrix_vec_product(m, &self.coeffs.to_vec())
    }
    pub fn mul_by_zq(&self, s: &Zq) -> Self {
        Self {
            param: self.param,
            coeffs: self.coeffs.iter().map(|c_i| *c_i * *s).collect(),
            evals: None,
        }
    }
    pub fn mul_by_u64(&self, s: u64) -> Self {
        let s = Zq::from_u64(self.param.q, s);
        Self {
            param: self.param,
            coeffs: self.coeffs.iter().map(|&e| e * s).collect(),
            evals: None,
        }
    }
    pub fn mul_by_f64(&self, s: f64) -> Self {
        Self {
            param: self.param,
            coeffs: self
                .coeffs
                .iter()
                .map(|c_i| Zq::from_f64(self.param.q, c_i.v as f64 * s))
                .collect(),
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
            .map(|e| (e.v as f64 / s as f64).round())
            .collect();
        Rq::from_vec_f64(&self.param, r)
    }

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // TODO simplify
        let mut str = "";
        let mut zero = true;
        for (i, coeff) in self.coeffs.iter().enumerate().rev() {
            if coeff.v == 0 {
                continue;
            }
            zero = false;
            f.write_str(str)?;
            if coeff.v != 1 {
                f.write_str(coeff.v.to_string().as_str())?;
                if i > 0 {
                    f.write_str("*")?;
                }
            }
            if coeff.v == 1 && i == 0 {
                f.write_str(coeff.v.to_string().as_str())?;
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
        f.write_str(self.param.q.to_string().as_str())?;
        f.write_str("/(X^")?;
        f.write_str(self.param.n.to_string().as_str())?;
        f.write_str("+1)")?;
        Ok(())
    }

    pub fn infinity_norm(&self) -> u64 {
        self.coeffs()
            .iter()
            .map(|x| {
                if x.v > (self.param.q / 2) {
                    self.param.q - x.v
                } else {
                    x.v
                }
            })
            .fold(0, |a, b| a.max(b))
    }
    pub fn mod_centered_q(&self) -> crate::ring_n::R {
        self.clone().to_r().mod_centered_q(self.param.q)
    }
}
pub fn matrix_vec_product(m: &Vec<Vec<Zq>>, v: &Vec<Zq>) -> Result<Vec<Zq>> {
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

    assert_eq!(m[0][0].q, v[0].q); // TODO change to returning err

    Ok(m.iter()
        .map(|row| {
            row.iter()
                .zip(v.iter())
                .map(|(&row_i, &v_i)| row_i * v_i)
                .sum()
        })
        .collect::<Vec<Zq>>())
}
pub fn transpose(m: &[Vec<Zq>]) -> Vec<Vec<Zq>> {
    assert!(m.len() > 0);
    assert!(m[0].len() > 0);
    let q = m[0][0].q;
    // TODO case when m[0].len()=0
    // TODO non square matrix
    let mut r: Vec<Vec<Zq>> = vec![vec![Zq::zero(q); m[0].len()]; m.len()];
    for (i, m_row) in m.iter().enumerate() {
        for (j, m_ij) in m_row.iter().enumerate() {
            r[j][i] = *m_ij;
        }
    }
    r
}

impl PartialEq for Rq {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs && self.param == other.param
    }
}
impl Add<Rq> for Rq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.param, rhs.param);
        Self {
            param: self.param,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l + r)
                .collect(),
            evals: None,
        }
    }
}
impl Add<&Rq> for &Rq {
    type Output = Rq;

    fn add(self, rhs: &Rq) -> Self::Output {
        assert_eq!(self.param, rhs.param);
        Rq {
            param: self.param,
            coeffs: zip_eq(self.coeffs.clone(), rhs.coeffs.clone())
                .map(|(l, r)| l + r)
                .collect(),
            evals: None,
        }
    }
}
impl AddAssign for Rq {
    fn add_assign(&mut self, rhs: Self) {
        debug_assert_eq!(self.param, rhs.param);
        for i in 0..self.param.n {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl Sum<Rq> for Rq {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let first = iter.next().unwrap();
        iter.fold(first, |acc, x| acc + x)
    }
}

impl Sub<Rq> for Rq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        assert_eq!(self.param, rhs.param);
        Self {
            param: self.param,
            coeffs: zip_eq(self.coeffs, rhs.coeffs)
                .map(|(l, r)| l - r)
                .collect(),
            evals: None,
        }
    }
}
impl Sub<&Rq> for &Rq {
    type Output = Rq;

    fn sub(self, rhs: &Rq) -> Self::Output {
        debug_assert_eq!(self.param, rhs.param);
        Rq {
            param: self.param,
            coeffs: zip_eq(self.coeffs.clone(), rhs.coeffs.clone())
                .map(|(l, r)| l - r)
                .collect(),
            evals: None,
        }
    }
}
impl SubAssign for Rq {
    fn sub_assign(&mut self, rhs: Self) {
        debug_assert_eq!(self.param, rhs.param);
        for i in 0..self.param.n {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl Mul<Rq> for Rq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        mul(&self, &rhs)
    }
}
impl Mul<&Rq> for &Rq {
    type Output = Rq;

    fn mul(self, rhs: &Rq) -> Self::Output {
        mul(self, rhs)
    }
}

// mul by Zq element
impl Mul<Zq> for Rq {
    type Output = Self;

    fn mul(self, s: Zq) -> Self {
        self.mul_by_zq(&s)
    }
}
impl Mul<&Zq> for &Rq {
    type Output = Rq;

    fn mul(self, s: &Zq) -> Self::Output {
        self.mul_by_zq(s)
    }
}
// mul by u64
impl Mul<u64> for Rq {
    type Output = Self;

    fn mul(self, s: u64) -> Self {
        self.mul_by_u64(s)
    }
}
impl Mul<&u64> for &Rq {
    type Output = Rq;

    fn mul(self, s: &u64) -> Self::Output {
        self.mul_by_u64(*s)
    }
}
// mul by f64
impl Mul<f64> for Rq {
    type Output = Self;

    fn mul(self, s: f64) -> Self {
        self.mul_by_f64(s)
    }
}
impl Mul<&f64> for &Rq {
    type Output = Rq;

    fn mul(self, s: &f64) -> Self::Output {
        self.mul_by_f64(*s)
    }
}

impl Neg for Rq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            param: self.param,
            coeffs: self.coeffs.iter().map(|c_i| -*c_i).collect(),
            evals: None,
        }
    }
}

// note: this assumes that Q is prime
fn mul_mut(lhs: &mut Rq, rhs: &mut Rq) -> Rq {
    assert_eq!(lhs.param, rhs.param);

    // reuse evaluations if already computed
    if !lhs.evals.is_some() {
        lhs.evals = Some(NTT::ntt(lhs).coeffs);
    };
    if !rhs.evals.is_some() {
        rhs.evals = Some(NTT::ntt(rhs).coeffs);
    };
    let lhs_evals = lhs.evals.clone().unwrap();
    let rhs_evals = rhs.evals.clone().unwrap();

    let c_ntt: Rq = Rq::from_vec(
        &lhs.param,
        zip_eq(lhs_evals, rhs_evals).map(|(l, r)| l * r).collect(),
    );
    let c = NTT::intt(&c_ntt);
    Rq::new(&lhs.param, c.coeffs, Some(c_ntt.coeffs))
}
// note: this assumes that Q is prime
// TODO impl karatsuba for non-prime Q. Alternatively check NTT with RNS trick.
fn mul(lhs: &Rq, rhs: &Rq) -> Rq {
    assert_eq!(lhs.param, rhs.param);

    // reuse evaluations if already computed
    let lhs_evals: Vec<Zq> = if lhs.evals.is_some() {
        lhs.evals.clone().unwrap()
    } else {
        NTT::ntt(lhs).coeffs
    };
    let rhs_evals: Vec<Zq> = if rhs.evals.is_some() {
        rhs.evals.clone().unwrap()
    } else {
        NTT::ntt(rhs).coeffs
    };

    let c_ntt: Rq = Rq::from_vec(
        &lhs.param,
        zip_eq(lhs_evals, rhs_evals).map(|(l, r)| l * r).collect(),
    );
    let c = NTT::intt(&c_ntt);
    Rq::new(&lhs.param, c.coeffs, Some(c_ntt.coeffs))
}

impl fmt::Display for Rq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)?;
        Ok(())
    }
}
impl fmt::Debug for Rq {
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
        let param = RingParam { q: 7, n: 3 };

        // p = 1x + 2x^2 + 3x^3 + 4 x^4 + 5 x^5 in R=Z_q[X]/(X^n +1)
        let p = Rq::from_vec_u64(&param, vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        // try with coefficients bigger than Q
        let p = Rq::from_vec_u64(&param, vec![0u64, 1, param.q + 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        // try with other ring
        let p = Rq::from_vec_u64(&RingParam { q: 7, n: 4 }, vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(p.to_string(), "3*x^3 + 2*x^2 + 3*x + 3 mod Z_7/(X^4+1)");

        let p = Rq::from_vec_u64(&param, vec![0u64, 0, 0, 0, 4, 5]);
        assert_eq!(p.to_string(), "2*x^2 + 3*x mod Z_7/(X^3+1)");

        let p = Rq::from_vec_u64(&param, vec![5u64, 4, 5, 2, 1, 0]);
        assert_eq!(p.to_string(), "5*x^2 + 3*x + 3 mod Z_7/(X^3+1)");

        let a = Rq::from_vec_u64(&param, vec![0u64, 1, 2, 3, 4, 5]);
        assert_eq!(a.to_string(), "4*x^2 + 4*x + 4 mod Z_7/(X^3+1)");

        let b = Rq::from_vec_u64(&param, vec![5u64, 4, 3, 2, 1, 0]);
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
        let param = RingParam {
            q: 2u64.pow(16) + 1,
            n: 4,
        };

        let a: Vec<u64> = vec![1u64, 2, 3, 4];
        let b: Vec<u64> = vec![1u64, 2, 3, 4];
        let c: Vec<u64> = vec![65513, 65517, 65531, 20];
        test_mul_opt(&param, a, b, c)?;

        let a: Vec<u64> = vec![0u64, 0, 0, 2];
        let b: Vec<u64> = vec![0u64, 0, 0, 2];
        let c: Vec<u64> = vec![0u64, 0, 65533, 0];
        test_mul_opt(&param, a, b, c)?;

        // TODO more testvectors

        Ok(())
    }
    fn test_mul_opt(
        param: &RingParam,
        a: Vec<u64>,
        b: Vec<u64>,
        expected_c: Vec<u64>,
    ) -> Result<()> {
        assert_eq!(a.len(), param.n);
        assert_eq!(b.len(), param.n);

        let mut a = Rq::from_vec_u64(&param, a);
        let mut b = Rq::from_vec_u64(&param, b);
        let expected_c = Rq::from_vec_u64(&param, expected_c);

        let c = mul_mut(&mut a, &mut b);
        assert_eq!(c, expected_c);
        Ok(())
    }

    #[test]
    fn test_rq_decompose() -> Result<()> {
        let param = RingParam { q: 16, n: 4 };
        let beta = 4;
        let l = 2;

        let a = Rq::from_vec_u64(&param, vec![7u64, 14, 3, 6]);
        let d = a.decompose(beta, l);

        assert_eq!(
            d[0].coeffs(),
            vec![1u64, 3, 0, 1]
                .iter()
                .map(|e| Zq::from_u64(param.q, *e))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            d[1].coeffs(),
            vec![3u64, 2, 3, 2]
                .iter()
                .map(|e| Zq::from_u64(param.q, *e))
                .collect::<Vec<_>>()
        );
        Ok(())
    }
}
