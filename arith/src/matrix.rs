use anyhow::{anyhow, Result};
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T>(pub Vec<Vec<T>>);

impl<T> Matrix<T>
where
    T: Copy + Add<Output = T> + Mul<Output = T> + Default + std::fmt::Debug,
{
    // TODO maybe rm this method, move it to tests only
    pub fn new(rows: usize, cols: usize, value: T) -> Self {
        Matrix(vec![vec![value; cols]; rows])
    }

    pub fn add(&self, other: &Matrix<T>) -> Result<Matrix<T>> {
        if self.0.len() != other.0.len() || self.0[0].len() != other.0[0].len() {
            return Err(anyhow!("dimensions don't match"));
        }

        let r = self
            .0
            .iter()
            .zip(&other.0)
            .map(|(row1, row2)| {
                row1.iter()
                    .zip(row2)
                    .map(|(a, b)| *a + *b)
                    .collect::<Vec<T>>()
            })
            .collect::<Vec<Vec<T>>>();

        Ok(Matrix(r))
    }

    pub fn mul(&self, other: &Matrix<T>) -> Result<Matrix<T>> {
        let rows_a = self.0.len();
        let cols_a = self.0[0].len();
        let rows_b = other.0.len();
        let cols_b = other.0[0].len();

        if cols_a != rows_b {
            return Err(anyhow!("self.n_cols != other.n_rows"));
        }

        let mut r = vec![vec![T::default(); cols_b]; rows_a];

        for i in 0..rows_a {
            for j in 0..cols_b {
                for k in 0..cols_a {
                    r[i][j] = r[i][j] + self.0[i][k] * other.0[k][j];
                }
            }
        }

        Ok(Matrix(r))
    }
    pub fn mul_vec(&self, v: &Vec<T>) -> Result<Vec<T>> {
        let rows = self.0.len();
        let cols = self.0[0].len();

        if cols != v.len() {
            return Err(anyhow!(
                "Number of columns in matrix does not match the length of the vector"
            ));
        }

        let mut r = vec![T::default(); rows];

        for i in 0..rows {
            for j in 0..cols {
                r[i] = r[i] + self.0[i][j] * v[j];
            }
        }
        Ok(r)
    }

    pub fn transpose(&self) -> Matrix<T> {
        let rows = self.0.len();
        let cols = self.0[0].len();
        let mut r = vec![vec![T::default(); rows]; cols];

        for i in 0..rows {
            for j in 0..cols {
                r[j][i] = self.0[i][j];
            }
        }

        Matrix(r)
    }

    pub fn scalar_mul(&self, scalar: T) -> Matrix<T> {
        let r = self
            .0
            .iter()
            .map(|row| row.iter().map(|&val| val * scalar).collect::<Vec<T>>())
            .collect::<Vec<Vec<T>>>();

        Matrix(r)
    }
}

// WIP. Currently uses ndarray, ndarray_linalg, num_complex to solve a system of
// linear equations A*x=b for x.
use crate::C;
impl Matrix<C<f64>> {
    pub fn solve(&self, b: &Vec<C<f64>>) -> Result<Vec<C<f64>>> {
        use ndarray::{Array1, Array2};
        use ndarray_linalg::Solve;
        use num_complex::Complex64;

        let m: Array2<Complex64> = Array2::from_shape_vec(
            (self.0.len(), self.0[0].len()),
            self.0
                .clone()
                .into_iter()
                .flatten()
                .map(|e| Complex64::new(e.re, e.im))
                .collect(),
        )
        .unwrap();
        let v: Array1<Complex64> = Array1::from_shape_vec(
            b.len(),
            b.iter().map(|e| Complex64::new(e.re, e.im)).collect(),
        )
        .unwrap();
        let r = m.solve(&v)?;
        let r: Vec<C<f64>> = r.into_iter().map(|e| C::<f64>::new(e.re, e.im)).collect();
        Ok(r)
    }
}
impl Matrix<f64> {
    // tmp (rm)
    pub fn solve(&self, b: Vec<f64>) -> Result<Vec<f64>> {
        use ndarray::{Array1, Array2};
        use ndarray_linalg::Solve;

        let m: Array2<f64> = Array2::from_shape_vec(
            (self.0.len(), self.0[0].len()),
            self.0.clone().into_iter().flatten().collect(),
        )
        .unwrap();
        let v: Array1<f64> = Array1::from_shape_vec(b.len(), b).unwrap();
        let r = m.solve(&v)?;
        let r: Vec<f64> = r.into_iter().map(|e| e).collect();
        Ok(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() -> Result<()> {
        let a = Matrix::new(2, 3, 1);
        let b = Matrix::new(2, 3, 2);
        let expected = Matrix::new(2, 3, 3);
        assert_eq!(a.add(&b).unwrap(), expected);
        Ok(())
    }

    #[test]
    fn test_mul() -> Result<()> {
        let a = Matrix::new(2, 3, 1);
        let b = Matrix::new(3, 2, 1);
        let expected = Matrix::new(2, 2, 3); // 2x3 * 3x2 = 2x2 matrix (with all values 3)
        assert_eq!(a.mul(&b).unwrap(), expected);
        Ok(())
    }

    #[test]
    fn test_transpose() -> Result<()> {
        let a = Matrix::new(2, 3, 1);
        let expected = Matrix::new(3, 2, 1);
        assert_eq!(a.transpose(), expected);
        Ok(())
    }

    #[test]
    fn test_scalar_mul() -> Result<()> {
        let a = Matrix::new(2, 3, 1);
        let expected = Matrix::new(2, 3, 3);
        assert_eq!(a.scalar_mul(3), expected);
        Ok(())
    }
}
