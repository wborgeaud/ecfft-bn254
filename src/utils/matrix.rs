use ark_ff::PrimeField;

#[derive(Clone, Copy)]
pub struct Matrix<T>(pub [[T; 2]; 2]);

impl<F: PrimeField> Matrix<F> {
    pub fn inverse(&self) -> Self {
        let [[a, b], [c, d]] = self.0;
        let det = a * d - b * c;
        Self([[d / det, -b / det], [-c / det, a / det]])
    }

    #[allow(clippy::many_single_char_names)]
    pub fn multiply(&self, v: [F; 2]) -> [F; 2] {
        let [[a, b], [c, d]] = self.0;
        let [x, y] = v;
        [a * x + b * y, c * x + d * y]
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{rand::Rng, test_rng};

    use crate::bn254::F;

    use super::Matrix;

    #[test]
    fn test_inverse() {
        let mut rng = test_rng();
        for _ in 0..100 {
            let a: F = rng.gen();
            let b: F = rng.gen();
            let c: F = rng.gen();
            let d: F = rng.gen();
            let mat = Matrix([[a, b], [c, d]]);
            let mat_inv = mat.inverse();
            let x: F = rng.gen();
            let y: F = rng.gen();
            let v = [x, y];

            assert_eq!(v, mat_inv.multiply(mat.multiply(v)));
            assert_eq!(v, mat.multiply(mat_inv.multiply(v)));
        }
    }
}
