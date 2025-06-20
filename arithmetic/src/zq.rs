use std::fmt;
use std::ops;

// Z_q, integers modulus q, not necessarily prime
#[derive(Clone, Copy, PartialEq)]
pub struct Zq<const Q: u64>(pub u64);

// WIP
// impl<const Q: u64> From<Vec<u64>> for Vec<Zq<Q>> {
//     fn from(v: Vec<u64>) -> Self {
//         v.into_iter().map(Zq::new).collect()
//     }
// }

impl<const Q: u64> Zq<Q> {
    pub fn new(e: u64) -> Self {
        if e >= Q {
            return Zq(e % Q);
        }
        Zq(e)
    }
    pub fn from_f64(e: f64) -> Self {
        // WIP method
        let e: i64 = e.round() as i64;
        if e < 0 {
            return Zq((Q as i64 + e) as u64);
        } else if e >= Q as i64 {
            return Zq((e % Q as i64) as u64);
        }
        Zq(e as u64)
    }
    pub fn from_bool(b: bool) -> Self {
        if b {
            Zq(1)
        } else {
            Zq(0)
        }
    }
    pub fn zero() -> Self {
        Zq(0u64)
    }
    pub fn square(self) -> Self {
        self * self
    }
    // modular exponentiation
    pub fn exp(self, e: Self) -> Self {
        // mul-square approach
        let mut res = Self(1);
        let mut rem = e.clone();
        let mut exp = self;
        // for rem != Self(0) {
        while rem != Self(0) {
            // if odd
            // TODO use a more readible expression
            if 1 - ((rem.0 & 1) << 1) as i64 == -1 {
                res = res * exp;
            }
            exp = exp.square();
            rem = Self(rem.0 >> 1);
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
        let mut r = Q;
        let mut new_t = 0;
        let mut new_r = self.0.clone();
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
        return Zq::new(t);
    }
    pub fn inv(self) -> Zq<Q> {
        let (g, x, _) = Self::egcd(self.0 as i128, Q as i128);
        if g != 1 {
            // None
            panic!("E");
        } else {
            let q = Q as i128;
            Zq(((x % q + q) % q) as u64) // TODO maybe just Zq::new(x)
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
}

impl<const Q: u64> Zq<Q> {
    fn r#mod(self) -> Self {
        if self.0 >= Q {
            return Zq(self.0 % Q);
        }
        self
    }
}

impl<const Q: u64> ops::Add<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut r = self.0 + rhs.0;
        if r >= Q {
            r -= Q;
        }
        Zq(r)
    }
}
impl<const Q: u64> ops::Add<&Zq<Q>> for &Zq<Q> {
    type Output = Zq<Q>;

    fn add(self, rhs: &Zq<Q>) -> Self::Output {
        let mut r = self.0 + rhs.0;
        if r >= Q {
            r -= Q;
        }
        Zq(r)
    }
}
impl<const Q: u64> ops::AddAssign<Zq<Q>> for Zq<Q> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}
impl<const Q: u64> std::iter::Sum for Zq<Q> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Zq(0), |acc, x| acc + x)
    }
}
impl<const Q: u64> ops::Sub<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Zq<Q> {
        if self.0 >= rhs.0 {
            Zq(self.0 - rhs.0)
        } else {
            Zq((Q + self.0) - rhs.0)
        }
    }
}
impl<const Q: u64> ops::Sub<&Zq<Q>> for &Zq<Q> {
    type Output = Zq<Q>;

    fn sub(self, rhs: &Zq<Q>) -> Self::Output {
        if self.0 >= rhs.0 {
            Zq(self.0 - rhs.0)
        } else {
            Zq((Q + self.0) - rhs.0)
        }
    }
}
impl<const Q: u64> ops::SubAssign<Zq<Q>> for Zq<Q> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}
impl<const Q: u64> ops::Neg for Zq<Q> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Zq(Q - self.0)
    }
}
impl<const Q: u64> ops::Mul<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Zq<Q> {
        // TODO non-naive way
        Zq(((self.0 as u128 * rhs.0 as u128) % Q as u128) as u64)
        // Zq((self.0 * rhs.0) % Q)
    }
}
impl<const Q: u64> ops::Div<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn div(self, rhs: Self) -> Zq<Q> {
        // TODO non-naive way
        // Zq((self.0 / rhs.0) % Q)
        self * rhs.inv()
    }
}

impl<const Q: u64> fmt::Display for Zq<Q> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl<const Q: u64> fmt::Debug for Zq<Q> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exp() {
        const Q: u64 = 1021;
        let a = Zq::<Q>(3);
        let b = Zq::<Q>(3);
        assert_eq!(a.exp(b), Zq(27));

        let a = Zq::<Q>(1000);
        let b = Zq::<Q>(3);
        assert_eq!(a.exp(b), Zq(949));
    }
    #[test]
    fn neg() {
        const Q: u64 = 1021;
        let a = Zq::<Q>::from_f64(101.0);
        let b = Zq::<Q>::from_f64(-1.0);
        assert_eq!(-a, a * b);
    }
}
