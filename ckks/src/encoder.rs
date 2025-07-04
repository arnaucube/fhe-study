use anyhow::Result;

use arith::{Matrix, Rq, C, R};

#[derive(Clone, Debug)]
pub struct SecretKey<const Q: u64, const N: usize>(Rq<Q, N>);

#[derive(Clone, Debug)]
pub struct PublicKey<const Q: u64, const N: usize>(Rq<Q, N>, Rq<Q, N>);

pub struct Encoder<const Q: u64, const N: usize> {
    scale_factor: C<f64>, // Δ (delta)
    primitive: C<f64>,
    basis: Matrix<C<f64>>,
    basis_t: Matrix<C<f64>>, // transposed basis
}

/// returns the mitive root of unity
fn primitive_root_of_unity(m: usize) -> C<f64> {
    let pi = C::<f64>::from(std::f64::consts::PI);
    ((C::<f64>::from(2f64) * pi * C::<f64>::i()) / C::<f64>::new(m as f64, 0f64)).exp()
}

/// where 'w' is 'omega', the primitive root of unity
fn vandermonde(n: usize, w: C<f64>) -> Matrix<C<f64>> {
    let mut v: Vec<Vec<C<f64>>> = vec![];
    for i in 0..n {
        let root = w.pow(2 * i as u32 + 1);
        let mut row: Vec<C<f64>> = vec![];
        for j in 0..n {
            row.push(root.pow(j as u32));
        }
        v.push(row);
    }
    Matrix::<C<f64>>(v)
}
impl<const Q: u64, const N: usize> Encoder<Q, N> {
    pub fn new(scale_factor: C<f64>) -> Self {
        let primitive: C<f64> = primitive_root_of_unity(2 * N);
        let basis = vandermonde(N, primitive);
        let basis_t = basis.transpose();

        Self {
            scale_factor,
            primitive,
            basis,
            basis_t,
        }
    }

    /// encode as described in the CKKS paper.
    /// from $\mathbb{C}^{N/2} \longrightarrow \mathbb{Z_q}[X]/(X^N +1) = R$
    // TODO use alg.1 from 2018-1043,
    //      or as in 2018-1073: $f(x) = 1N (U^T.conj() m + U^T m.conj())$
    pub fn encode(&self, z: &[C<f64>]) -> Result<R<N>> {
        // $pi^{-1}: \mathbb{C}^{N/2} \longrightarrow \mathbb{H}$
        let expanded = self.pi_inv(z);

        // scale the values
        let scaled: Vec<C<f64>> = expanded.iter().map(|e| *e * self.scale_factor).collect();

        // but $\mathbb{H} \neq \sigma(R)$, since $\sigma(R) \subseteq \mathbb{H}$, so we need to
        // discretize $\pi^{-1}(z)$ into an element of $\sigma(R)$.

        // discretize \pi^-1(z_projected) to \sigma(R)
        // project 'scaled' into \sigma(R):
        // get the orthogonal basis (note: that would be doing Gram-Schmidt, which is not this, but
        // we're fine since the basis=Vandermonde matrix which is orthogonal, so we project z to it):
        // $z = \sum z_i * b_i, with z_i = <z,b_i>/||b_i||^2$
        let z_projected = self
            .basis_t
            .0
            .iter()
            .map(|b_i| {
                // TODO: the b_j.conj() can be precomputed at initialization (of the basis)
                let num: C<f64> = scaled
                    .iter()
                    .zip(b_i.iter())
                    .map(|(z_j, b_j)| *z_j * b_j.conj())
                    .sum::<C<f64>>();
                let den: C<f64> = b_i.iter().map(|b_j| *b_j * b_j.conj()).sum::<C<f64>>();
                let mut z_i = num / den;
                z_i.im = 0.0; // get only the real component
                z_i
            })
            .collect::<Vec<C<f64>>>();

        // V * z_projected (V: Vandermonde matrix)
        let discretized = self.basis.mul_vec(&z_projected)?;

        // sigma_inv
        let r = self.sigma_inv(&discretized)?;

        // TMP: naive round, maybe do gaussian
        let coeffs = r.iter().map(|e| e.re.round() as i64).collect::<Vec<i64>>();
        Ok(R::from_vec(coeffs))
    }

    pub fn decode(&self, p: &R<N>) -> Result<Vec<C<f64>>> {
        let p: Vec<C<f64>> = p
            .coeffs()
            .iter()
            .map(|&e| C::<f64>::new(e as f64, 0_f64)) // TODO review u64 to f64 conversion overflow
            .collect();
        let in_sigma = self.sigma(&p)?;

        let deescalated: Vec<C<f64>> = in_sigma.iter().map(|e| *e / self.scale_factor).collect();
        Ok(self.pi(&deescalated))
    }

    /// pi: \mathbb{H} \longrightarrow \mathbb{C}^{N/2}
    fn pi(&self, z: &[C<f64>]) -> Vec<C<f64>> {
        z[..N / 2].to_vec()
    }
    /// pi^{-1}: \mathbb{C}^{N/2} \longrightarrow \mathbb{H}
    fn pi_inv(&self, z: &[C<f64>]) -> Vec<C<f64>> {
        z.iter()
            .cloned()
            .chain(z.iter().rev().map(|z_i| z_i.conj()))
            .collect()
    }

    fn sigma(&self, p: &[C<f64>]) -> Result<Vec<C<f64>>> {
        // the roots of unity are already calculated in the 2nd row of the transpose of the
        // Vandermonde matrix used as the basis (ie. the 2nd column of the Vandermonde matrix).
        // let roots = &self.basis_t[1];
        // // Approach 1: evaluate p at the roots of unity
        // let mut z = vec![];
        // for root_i in roots.iter() {
        //     z.push(eval(p, root_i));
        // }

        // Approach 2: Vandermonde * p
        let z: Vec<C<f64>> = self.basis.mul_vec(&p.to_vec())?;

        // TODO check using NTT-ish (2018-1043) for the encode/decode

        Ok(z)
    }

    fn sigma_inv(&self, z: &Vec<C<f64>>) -> Result<Vec<C<f64>>> {
        // $\alpha = A^{-1} * z$
        let a = self.basis.solve(z)?;
        Ok(a.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_encode_decode() -> Result<()> {
        const Q: u64 = 1024;
        // const N: usize = 4; // ie. m=2*n=8
        const N: usize = 16;

        let T = 16; // WIP
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let z: Vec<C<f64>> = std::iter::repeat_with(|| {
                C::<f64>::new(rng.gen_range(0..T) as f64, rng.gen_range(0..T) as f64)
            })
            .take(N / 2)
            .collect();

            let delta = C::<f64>::new(64.0, 0.0); // delta = scaling factor
            let encoder = Encoder::<Q, N>::new(delta);

            let m: R<N> = encoder.encode(&z)?; // polynomial (encoded vec) \in R

            let z_decoded = encoder.decode(&m)?;

            // round it to compare it to the initial value
            let rounded_z_decoded: Vec<C<f64>> = z_decoded
                .iter()
                .map(|c| C::<f64>::new(c.re.round(), c.im.round()))
                .collect();
            assert_eq!(rounded_z_decoded, z);
        }

        Ok(())
    }
}
