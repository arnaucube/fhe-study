use rand::{distributions::Distribution, Rng};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Z_q, integers modulus q, not necessarily prime
#[derive(Clone, Copy, PartialEq)]
pub struct Zq<const Q: u64>(pub u64);

// WIP
// impl<const Q: u64> From<Vec<u64>> for Vec<Zq<Q>> {
//     fn from(v: Vec<u64>) -> Self {
//         v.into_iter().map(Zq::new).collect()
//     }
// }

pub(crate) fn modulus_u64<const Q: u64>(e: u64) -> u64 {
    (e % Q + Q) % Q
}
impl<const Q: u64> Zq<Q> {
    pub fn rand(mut rng: impl Rng, dist: impl Distribution<f64>) -> Self {
        // TODO WIP
        let r: f64 = dist.sample(&mut rng);
        Self::from_f64(r)
        // Self::from_u64(r.round() as u64)
    }
    pub fn from_u64(e: u64) -> Self {
        if e >= Q {
            // (e % Q + Q) % Q
            return Zq(modulus_u64::<Q>(e));
            // return Zq(e % Q);
        }
        Zq(e)
    }
    pub fn from_f64(e: f64) -> Self {
        // WIP method
        let e: i64 = e.round() as i64;
        let q = Q as i64;
        if e < 0 || e >= q {
            return Zq(((e % q + q) % q) as u64);
        }
        Zq(e as u64)

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
    pub fn from_bool(b: bool) -> Self {
        if b {
            Zq(1)
        } else {
            Zq(0)
        }
    }
    pub fn zero() -> Self {
        Self(0u64)
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
        return Zq::from_u64(t);
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

    /// perform the mod switch operation from Q to Q', where Q2=Q'
    pub fn mod_switch<const Q2: u64>(&self) -> Zq<Q2> {
        Zq::<Q2>::from_u64(((self.0 as f64 * Q2 as f64) / Q as f64).round() as u64)
    }

    pub fn decompose(&self, beta: u32, l: u32) -> Vec<Self> {
        if beta == 2 {
            self.decompose_base2(l)
        } else {
            self.decompose_base_beta(beta, l)
        }
    }
    pub fn decompose_base_beta(&self, beta: u32, l: u32) -> Vec<Self> {
        let mut rem: u64 = self.0;
        // next if is for cases in which beta does not divide Q (concretely
        // beta^l!=Q). round to the nearest multiple of q/beta^l
        if rem >= beta.pow(l) as u64 {
            // rem = Q - 1 - (Q / beta as u64); // floor
            return vec![Zq(beta as u64 - 1); l as usize];
        }

        let mut x: Vec<Self> = vec![];
        for i in 1..l + 1 {
            let den = Q / beta.pow(i) as u64;
            let x_i = rem / den; // division between u64 already does floor
            x.push(Self::from_u64(x_i));
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
        if self.0 >= 1 << l as u64 {
            // rem = Q - 1 - (Q / beta as u64); // floor
            // (where beta=2)
            return vec![Zq(1); l as usize];
        }

        (0..l)
            .rev()
            .map(|i| Self(((self.0 >> i) & 1) as u64))
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

impl<const Q: u64> Zq<Q> {
    fn r#mod(self) -> Self {
        if self.0 >= Q {
            return Zq(self.0 % Q);
        }
        self
    }
}

impl<const Q: u64> Add<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut r = self.0 + rhs.0;
        if r >= Q {
            r -= Q;
        }
        Zq(r)
    }
}
impl<const Q: u64> Add<&Zq<Q>> for &Zq<Q> {
    type Output = Zq<Q>;

    fn add(self, rhs: &Zq<Q>) -> Self::Output {
        let mut r = self.0 + rhs.0;
        if r >= Q {
            r -= Q;
        }
        Zq(r)
    }
}
impl<const Q: u64> AddAssign<Zq<Q>> for Zq<Q> {
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
impl<const Q: u64> Sub<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Zq<Q> {
        if self.0 >= rhs.0 {
            Zq(self.0 - rhs.0)
        } else {
            Zq((Q + self.0) - rhs.0)
        }
    }
}
impl<const Q: u64> Sub<&Zq<Q>> for &Zq<Q> {
    type Output = Zq<Q>;

    fn sub(self, rhs: &Zq<Q>) -> Self::Output {
        if self.0 >= rhs.0 {
            Zq(self.0 - rhs.0)
        } else {
            Zq((Q + self.0) - rhs.0)
        }
    }
}
impl<const Q: u64> SubAssign<Zq<Q>> for Zq<Q> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}
impl<const Q: u64> Neg for Zq<Q> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.0 == 0 {
            return self;
        }
        Zq(Q - self.0)
    }
}
impl<const Q: u64> Mul<Zq<Q>> for Zq<Q> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Zq<Q> {
        // TODO non-naive way
        Zq(((self.0 as u128 * rhs.0 as u128) % Q as u128) as u64)
        // Zq((self.0 * rhs.0) % Q)
    }
}
impl<const Q: u64> Div<Zq<Q>> for Zq<Q> {
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
    use rand::distributions::Uniform;

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

    fn recompose<const Q: u64>(beta: u32, l: u32, d: Vec<Zq<Q>>) -> Zq<Q> {
        let mut x = 0u64;
        for i in 0..l {
            x += d[i as usize].0 * Q / beta.pow(i + 1) as u64;
        }
        Zq::from_u64(x)
    }

    #[test]
    fn test_decompose() {
        const Q1: u64 = 16;
        let beta: u32 = 2;
        let l: u32 = 4;
        let x = Zq::<Q1>::from_u64(9);
        let d = x.decompose(beta, l);

        assert_eq!(recompose::<Q1>(beta, l, d), x);

        const Q: u64 = 5u64.pow(3);
        let beta: u32 = 5;
        let l: u32 = 3;

        let dist = Uniform::new(0_u64, Q);
        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let x = Zq::<Q>::from_u64(dist.sample(&mut rng));
            let d = x.decompose(beta, l);
            assert_eq!(d.len(), l as usize);
            assert_eq!(recompose::<Q>(beta, l, d), x);
        }
    }

    #[test]
    fn test_decompose_approx() {
        const Q: u64 = 2u64.pow(4) + 1;
        let beta: u32 = 2;
        let l: u32 = 4;
        let x = Zq::<Q>::from_u64(16); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(recompose::<Q>(beta, l, d), Zq(15));

        const Q2: u64 = 5u64.pow(3) + 1;
        let beta: u32 = 5;
        let l: u32 = 3;
        let x = Zq::<Q2>::from_u64(125); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(recompose::<Q2>(beta, l, d), Zq(124));

        const Q3: u64 = 2u64.pow(16) + 1;
        let beta: u32 = 2;
        let l: u32 = 16;
        let x = Zq::<Q3>::from_u64(Q3 - 1); // in q, but bigger than beta^l
        let d = x.decompose(beta, l);
        assert_eq!(d.len(), l as usize);
        assert_eq!(recompose::<Q3>(beta, l, d), Zq(beta.pow(l) as u64 - 1));
    }
}
