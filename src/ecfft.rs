use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial};

use crate::utils::{isogeny::Isogeny, matrix::Matrix};

pub trait EcFftParameters<F: PrimeField>: Sized {
    const LOG_N: usize;
    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F>;

    fn isogenies() -> Vec<Isogeny<F>>;

    fn precompute() -> EcFftPrecomputation<F, Self>;
}

pub struct EcFftPrecomputationStep<F: PrimeField, P: EcFftParameters<F>> {
    pub s: Vec<F>,
    pub s_prime: Vec<F>,
    pub matrices: Vec<Matrix<F>>,
    pub inverse_matrices: Vec<Matrix<F>>,
    pub _phantom: PhantomData<P>,
}
pub struct EcFftPrecomputation<F: PrimeField, P: EcFftParameters<F>> {
    pub steps: Vec<EcFftPrecomputationStep<F, P>>,
    pub final_s: F,
    pub final_s_prime: F,
    pub _phantom: PhantomData<P>,
}

impl<F: PrimeField, P: EcFftParameters<F>> EcFftPrecomputation<F, P> {
    pub fn extend(&self, evals: &[F]) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }
        debug_assert_eq!(
            n & (n - 1),
            0,
            "The number of evaluations should be a power of 2."
        );
        let log_n = n.trailing_zeros() as usize;
        debug_assert!(
            log_n < P::LOG_N,
            "Got {} evaluations, can extend at most {} evaluations.",
            n,
            1 << (P::LOG_N - 1)
        );

        let EcFftPrecomputationStep {
            matrices,
            inverse_matrices,
            _phantom,
            ..
        } = &self.steps[P::LOG_N - 1 - log_n];
        let mut p0_evals = Vec::new();
        let mut p1_evals = Vec::new();
        let nn = n / 2;
        for j in 0..nn {
            let (y0, y1) = (evals[j], evals[j + nn]);
            let [p0, p1] = inverse_matrices[j].multiply([y0, y1]);
            p0_evals.push(p0);
            p1_evals.push(p1);
        }
        let p0_evals_prime = self.extend(&p0_evals);
        let p1_evals_prime = self.extend(&p1_evals);

        let mut ans_left = Vec::new();
        let mut ans_right = Vec::new();
        for (m, (&p0, &p1)) in matrices
            .iter()
            .zip(p0_evals_prime.iter().zip(&p1_evals_prime))
        {
            let [x, y] = m.multiply([p0, p1]);
            ans_left.push(x);
            ans_right.push(y);
        }

        let mut ans = ans_left;
        ans.append(&mut ans_right);
        ans
    }

    // L = S u S' -> s u s'
    //
    pub fn evaluate_over_domain(&self, poly: &DensePolynomial<F>, log_n: usize) -> Vec<F> {
        if log_n == 0 {
            return vec![poly.coeffs[0]];
        }
        let n = 1 << log_n;
        assert!(poly.len() <= n);
        let mut coeffs = poly.coeffs.clone();
        coeffs.resize(n, F::zero());
        let low = coeffs[..n / 2].to_vec();
        let high = coeffs[n / 2..].to_vec();
        let low_evals = self.evaluate_over_domain(&DensePolynomial { coeffs: low }, log_n - 1);
        let high_evals = self.evaluate_over_domain(&DensePolynomial { coeffs: high }, log_n - 1);
        let low_evals_prime = self.extend(&low_evals);
        let high_evals_prime = self.extend(&high_evals);
        dbg!(poly);
        dbg!(&low_evals, &high_evals, &low_evals_prime, &high_evals_prime);

        let EcFftPrecomputationStep { s, s_prime, .. } = &self.steps[P::LOG_N - 1 - log_n];
        debug_assert_eq!(s.len(), n);
        let mut ans = Vec::new();
        for i in 0..n / 2 {
            ans.push(low_evals[i] + s[2 * i].pow([n as u64 / 2]) * high_evals[i]);
            ans.push(low_evals_prime[i] + s[2 * i + 1].pow([n as u64 / 2]) * high_evals_prime[i]);
        }

        ans
    }
}

#[cfg(test)]
mod tests {
    use crate::bn254::{Bn254EcFftParameters, F};
    use ark_poly::{univariate::DensePolynomial, Polynomial};
    use ark_std::{
        rand::{distributions::Standard, prelude::Distribution, Rng},
        test_rng,
    };

    use super::*;

    fn test_extend_i<F: PrimeField, P: EcFftParameters<F>>(
        i: usize,
        precomputation: &EcFftPrecomputation<F, P>,
    ) where
        Standard: Distribution<F>,
    {
        let n = 1 << i;
        let mut rng = test_rng();
        let coeffs: Vec<F> = (0..n).map(|_| rng.gen()).collect();
        let poly = DensePolynomial { coeffs };
        let EcFftPrecomputationStep { s, s_prime, .. } =
            &precomputation.steps[Bn254EcFftParameters::LOG_N - 1 - i];
        let evals_s = s.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
        let evals_s_prime = s_prime.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
        assert_eq!(evals_s_prime, precomputation.extend(&evals_s));
    }

    #[test]
    fn test_extend() {
        let precomputation = Bn254EcFftParameters::precompute();
        for i in 1..Bn254EcFftParameters::LOG_N {
            test_extend_i::<F, _>(i, &precomputation);
        }
    }

    fn test_eval_i<F: PrimeField, P: EcFftParameters<F>>(
        i: usize,
        precomputation: &EcFftPrecomputation<F, P>,
    ) where
        Standard: Distribution<F>,
    {
        let n = 1 << i;
        let mut rng = test_rng();
        let coeffs: Vec<F> = (0..n).map(|_| rng.gen()).collect();
        let poly = DensePolynomial { coeffs };
        let EcFftPrecomputationStep { s, s_prime, .. } =
            &precomputation.steps[Bn254EcFftParameters::LOG_N - 1 - i];
        let domain = (0..n / 2)
            .map(|i| [s[i], s_prime[i]])
            .flatten()
            .collect::<Vec<_>>();
        let evals = s.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
        assert_eq!(evals, precomputation.evaluate_over_domain(&poly, i));
    }

    #[test]
    fn test_eval() {
        let precomputation = Bn254EcFftParameters::precompute();
        for i in 2..=Bn254EcFftParameters::LOG_N {
            test_eval_i::<F, _>(i, &precomputation);
        }
    }
}
