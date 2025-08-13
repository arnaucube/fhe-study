use rand::{distributions::Distribution, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Z_q, integers modulus q, not necessarily prime
#[derive(Clone, Copy, PartialEq)]
pub struct Zq {
    pub q: u64,
    pub v: u64,
}

// WIP
// impl<const Q: u64> From<Vec<u64>> for Vec<Zq<Q>> {
//     fn from(v: Vec<u64>) -> Self {
//         v.into_iter().map(Zq::new).collect()
//     }
// }

pub(crate) fn modulus_u64(q: u64, e: u64) -> u64 {
    (e % q + q) % q
}
impl Zq {
    pub fn rand(mut rng: impl Rng, dist: impl Distribution<f64>, q: u64) -> Self {
        // TODO WIP
        let r: f64 = dist.sample(&mut rng);
        Self::from_f64(q, r)
        // Self::from_u64(r.round() as u64)
    }
    pub fn from_u64(q: u64, v: u64) -> Self {
        if v >= q {
            // (v % Q + Q) % Q
            return Zq {
                q,
                v: modulus_u64(q, v),
            };
            // return Zq(v % Q);
        }
        Zq { q, v }
    }
    pub fn from_f64(q: u64, e: f64) -> Self {
        // WIP method
        let e: i64 = e.round() as i64;
        let q_i64 = q as i64;
        if e < 0 || e >= q_i64 {
            return Zq::from_u64(q, ((e % q_i64 + q_i64) % q_i64) as u64);
        }
        Zq { q, v: e as u64 }

        // if e < 0 {
        //     // dbg!(&e);
        //     // dbg!(Zq::<Q>(((Q as i64 + e) % Q as i64) as u64));
        //     // return Zq(((Q as i64 + e) % Q as i64) as u64);
        //     return Zq(e as u64 % Q);
        // } else if e >= Q as i64 {
        //     return Zq((e % Q as i64) as u64);
        // }
        // Zq(e as u64)
    }
    pub fn from_bool(q: u64, b: bool) -> Self {
        if b {
            Zq { q, v: 1 }
        } else {
            Zq { q, v: 0 }
        }
    }
    pub fn zero(q: u64) -> Self {
        Self { q, v: 0u64 }
    }
    pub fn one(q: u64) -> Self {
        Self { q, v: 1u64 }
    }
    pub fn square(self) -> Self {
        self * self
    }
    // modular exponentiation
    pub fn exp(self, e: Self) -> Self {
        // mul-square approach
        let mut res = Self::one(self.q);
        let mut rem = e.clone();
        let mut exp = self;
        // for rem != Self(0) {
        while rem != Self::zero(self.q) {
            // if odd
            // TODO use a more readible expression
            if 1 - ((rem.v & 1) << 1) as i64 == -1 {
                res = res * exp;
            }
            exp = exp.square();
            rem = Self {
                q: self.q,
                v: rem.v >> 1,
            };
        }
        res
    }
    // multiplicative inverse
    // WARNING: if this is needed, it means that 'Zq' is a Finite Field. For the moment we assume
    // we work in a Finite Field
    pub fn inv_OLD(self) -> Self {
        // TODO
        // let a = self.0;
        // let q = Q;
        let mut t = 0;
        let mut r = self.q;
        let mut new_t = 0;
        let mut new_r = self.v.clone();
        while new_r != 0 {
            let q = r / new_r;

            t = new_t.clone();
            new_t = t - q;

            r = new_r.clone();
            new_r = r - (q * new_r);
        }
        // if t < 0 {
        //     t = t + q;
        // }
        return Zq::from_u64(self.q, t);
    }
    pub fn inv(self) -> Zq {
        let (g, x, _) = Self::egcd(self.v as i128, self.q as i128);
        if g != 1 {
            // None
            panic!("E");
        } else {
            let q = self.q as i128;
            Zq::from_u64(self.q, ((x % q + q) % q) as u64) // TODO maybe just Zq::new(x)
        }
    }
    fn egcd(a: i128, b: i128) -> (i128, i128, i128) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (g, x, y) = Self::egcd(b % a, a);
            (g, y - (b / a) * x, x)
        }
    }

    /// perform the mod switch operation from Q to Q', where Q2=Q'
    pub fn mod_switch(&self, q2: u64) -> Zq {
        Zq::from_u64(
            q2,
            ((self.v as f64 * q2 as f64) / self.q as f64).round() as u64,
        )
    }

    pub fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        if beta == 2 {
            self.decompose_base2(l)
        } else {
            self.decompose_base_beta(beta, l)
        }
    }
    pub fn decompose_base_beta(&self, beta: u32, l: u32) -> Vec<Self> {
        let mut rem: u64 = self.v;
        // next if is for cases in which beta does not divide Q (concretely
        // beta^l!=Q). round to the nearest multiple of q/beta^l
        if rem >= beta.pow(l) as u64 {
            // rem = Q - 1 - (Q / beta as u64); // floor
            return vec![
                Zq {
                    q: self.q,
                    v: beta as u64 - 1
                };
                l as usize
            ];
        }

        let mut x: Vec<Self> = vec![];
        for i in 1..l + 1 {
            let den = self.q / beta.pow(i) as u64;
            let x_i = rem / den; // division between u64 already does floor
            x.push(Self::from_u64(self.q, x_i));
            if x_i != 0 {
                rem = rem % den;
            }
        }
        x
    }
    /// decompose when beta=2
    pub fn decompose_base2(&self, l: u32) -> Vec<Self> {
        // next if is for cases in which beta does not divide Q (concretely
        // beta^l!=Q). round to the nearest multiple of q/beta^l
        if self.v >= 1 << l as u64 {
            // rem = Q - 1 - (Q / beta as u64); // floor
            // (where beta=2)
            return vec![Zq::one(self.q); l as usize];
        }

        (0..l)
            .rev()
            .map(|i| Self::from_u64(self.q, ((self.v >> i) & 1) as u64))
            .collect()

        // naive ver:
        // let mut rem: u64 = self.0;
        // // next if is for cases in which beta does not divide Q (concretely
        // // beta^l!=Q). round to the nearest multiple of q/beta^l
        // if rem >= 1 << l as u64 {
        //     // rem = Q - 1 - (Q / beta as u64); // floor
        //     return vec![Zq(1); l as usize];
        // }
        //
        // let mut x: Vec<Self> = vec![];
        // for i in 1..l + 1 {
        //     let den = Q / (1 << i as u64);
        //     let x_i = rem / den; // division between u64 already does floor
        //     x.push(Self::from_u64(x_i));
        //     if x_i != 0 {
        //         rem = rem % den;
        //     }
        // }
        // x
    }
}

impl Zq {
    fn r#mod(self) -> Self {
        if self.v >= self.q {
            return Zq::from_u64(self.q, self.v % self.q);
        }
        self
    }
}

impl Add<Zq> for Zq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.q, rhs.q);

        let mut v = self.v + rhs.v;
        if v >= self.q {
            v -= self.q;
        }
        Zq { q: self.q, v }
    }
}
impl Add<&Zq> for &Zq {
    type Output = Zq;

    fn add(self, rhs: &Zq) -> Self::Output {
        assert_eq!(self.q, rhs.q);

        let mut v = self.v + rhs.v;
        if v >= self.q {
            v -= self.q;
        }
        Zq { q: self.q, v }
    }
}
impl AddAssign<Zq> for Zq {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}
impl std::iter::Sum for Zq {
    fn sum<I>(mut iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let first: Zq = iter.next().unwrap();
        iter.fold(first, |acc, x| acc + x)
    }
}
impl Sub<Zq> for Zq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Zq {
        assert_eq!(self.q, rhs.q);

        if self.v >= rhs.v {
            Zq {
                q: self.q,
                v: self.v - rhs.v,
            }
        } else {
            Zq {
                q: self.q,
                v: (self.q + self.v) - rhs.v,
            }
        }
    }
}
impl Sub<&Zq> for &Zq {
    type Output = Zq;

    fn sub(self, rhs: &Zq) -> Self::Output {
        assert_eq!(self.q, rhs.q);

        if self.q >= rhs.q {
            Zq {
                q: self.q,
                v: self.v - rhs.v,
            }
        } else {
            Zq {
                q: self.q,
                v: (self.q + self.v) - rhs.v,
            }
        }
    }
}
impl SubAssign<Zq> for Zq {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}
impl Neg for Zq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.v == 0 {
            return self;
        }
        Zq {
            q: self.q,
            v: self.q - self.v,
        }
    }
}
impl Mul<Zq> for Zq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Zq {
        assert_eq!(self.q, rhs.q);

        // TODO non-naive way
        Zq {
            q: self.q,
            v: ((self.v as u128 * rhs.v as u128) % self.q as u128) as u64,
        }
        // Zq((self.0 * rhs.0) % Q)
    }
}
impl Div<Zq> for Zq {
    type Output = Self;

    fn div(self, rhs: Self) -> Zq {
        // TODO non-naive way
        // Zq((self.0 / rhs.0) % Q)
        self * rhs.inv()
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.v)
    }
}
impl fmt::Debug for Zq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::Uniform;

    #[test]
    fn exp() {
        const q: u64 = 1021;
        let a = Zq::from_u64(q, 3);
        let b = Zq::from_u64(q, 3);
        assert_eq!(a.exp(b), Zq::from_u64(q, 27));

        let a = Zq::from_u64(q, 1000);
        let b = Zq::from_u64(q, 3);
        assert_eq!(a.exp(b), Zq::from_u64(q, 949));
    }
    #[test]
    fn neg() {
        let q: u64 = 1021;
        let a = Zq::from_f64(q, 101.0);
        let b = Zq::from_f64(q, -1.0);
        assert_eq!(-a, a * b);
    }

    fn recompose(q: u64, beta: u32, l: u32, d: Vec<Zq>) -> Zq {
        let mut x = 0u64;
        for i in 0..l {
            x += d[i as usize].v * q / beta.pow(i + 1) as u64;
        }
        Zq::from_u64(q, x)
    }

    #[test]
    fn test_decompose() {
        let q1: u64 = 16;
        let beta: u32 = 2;
        let l: u32 = 4;
        let x = Zq::from_u64(q1, 9);
        let d = x.decompose(beta, l);

        assert_eq!(recompose(q1, beta, l, d), x);

        let q: u64 = 5u64.pow(3);
        let beta: u32 = 5;
        let l: u32 = 3;

        let dist = Uniform::new(0_u64, q);
        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let x = Zq::from_u64(q, dist.sample(&mut rng));
            let d = x.decompose(beta, l);
            assert_eq!(d.len(), l as usize);
            assert_eq!(recompose(q, beta, l, d), x);
        }
    }

    #[test]
    fn test_decompose_approx() {
        let q: u64 = 2u64.pow(4) + 1;
        let beta: u32 = 2;
        let l: u32 = 4;
        let x = Zq::from_u64(q, 16); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(recompose(q, beta, l, d), Zq::from_u64(q, 15));

        let q2: u64 = 5u64.pow(3) + 1;
        let beta: u32 = 5;
        let l: u32 = 3;
        let x = Zq::from_u64(q2, 125); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(recompose(q2, beta, l, d), Zq::from_u64(q2, 124));

        let q3: u64 = 2u64.pow(16) + 1;
        let beta: u32 = 2;
        let l: u32 = 16;
        let x = Zq::from_u64(q3, q3 - 1); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(
            recompose(q3, beta, l, d),
            Zq::from_u64(q3, beta.pow(l) as u64 - 1)
        );
    }
}
