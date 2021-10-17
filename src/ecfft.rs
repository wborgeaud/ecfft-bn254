use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;

use crate::utils::{isogeny::Isogeny, matrix::Matrix};

pub trait EcFftParameters<F: PrimeField>: Sized {
    /// Logarithm of the size of the maximal ECFFT coset in the curve.
    const LOG_N: usize;
    /// Size of the maximal ECFFT coset in the curve.
    const N: usize = 1 << Self::LOG_N;

    /// Maximal ECFFT coset in the curve.
    fn coset() -> Vec<F>;

    /// `Self::coset()[::i]`
    fn sub_coset(i: usize) -> Vec<F> {
        Self::coset().into_iter().step_by(1 << i).collect()
    }

    /// Isogenies of degree 2 used to compute the ECFFT.
    /// They form a chain of isogenies `E -> E' -> E'' -> ...` starting from the orignal curve.
    fn isogenies() -> Vec<Isogeny<F>>;

    /// Computes the ECFFT precomputations on a given coset.
    fn precompute_on_coset(coset: &[F]) -> EcFftCosetPrecomputation<F, Self> {
        let n = coset.len();
        let log_n = n.trailing_zeros() as usize;
        let isogenies = Self::isogenies();
        debug_assert_eq!(isogenies.len(), Self::LOG_N - 1);

        let mut s = coset.iter().step_by(2).copied().collect::<Vec<_>>();
        let mut s_prime = coset.iter().skip(1).step_by(2).copied().collect::<Vec<_>>();

        let mut steps = Vec::new();
        for i in (1..log_n).rev() {
            let n = 1 << i;
            let nn = n / 2;
            let q = nn - 1;
            let psi = isogenies[log_n - 1 - i];
            let mut matrices = Vec::new();
            let mut inverse_matrices = Vec::new();
            for j in 0..nn {
                let (s0, s1) = (s[j], s[j + nn]);
                debug_assert_eq!(psi.eval(s0), psi.eval(s1));
                inverse_matrices.push(
                    Matrix([
                        [
                            psi.eval_den(s0).pow([q as u64]),
                            s0 * psi.eval_den(s0).pow([q as u64]),
                        ],
                        [
                            psi.eval_den(s1).pow([q as u64]),
                            s1 * psi.eval_den(s1).pow([q as u64]),
                        ],
                    ])
                    .inverse(),
                );

                let (s0, s1) = (s_prime[j], s_prime[j + nn]);
                debug_assert_eq!(psi.eval(s0), psi.eval(s1));
                matrices.push(Matrix([
                    [
                        psi.eval_den(s0).pow([q as u64]),
                        s0 * psi.eval_den(s0).pow([q as u64]),
                    ],
                    [
                        psi.eval_den(s1).pow([q as u64]),
                        s1 * psi.eval_den(s1).pow([q as u64]),
                    ],
                ]));
            }
            steps.push(EcFftPrecomputationStep::<F, Self> {
                s: s.clone(),
                s_prime: s_prime.clone(),
                matrices,
                inverse_matrices,
                _phantom: PhantomData,
            });
            s = s.into_iter().take(nn).map(|x| psi.eval(x)).collect();
            s_prime = s_prime.into_iter().take(nn).map(|x| psi.eval(x)).collect();
        }
        debug_assert_eq!((s.len(), s_prime.len()), (1, 1));

        EcFftCosetPrecomputation {
            coset: coset.to_vec(),
            steps,
        }
    }

    /// Computes the ECFFT precomputations of all `Self::sub_coset(i)`.
    fn precompute() -> EcFftPrecomputation<F, Self> {
        let mut coset_precomputations = Vec::new();
        let mut coset = Self::coset();
        for _ in 0..Self::LOG_N {
            coset_precomputations.push(Self::precompute_on_coset(&coset));
            coset = coset.into_iter().step_by(2).collect();
        }
        EcFftPrecomputation {
            coset_precomputations,
        }
    }
}

pub struct EcFftPrecomputationStep<F: PrimeField, P: EcFftParameters<F>> {
    pub s: Vec<F>,
    pub s_prime: Vec<F>,
    pub matrices: Vec<Matrix<F>>,
    pub inverse_matrices: Vec<Matrix<F>>,
    pub _phantom: PhantomData<P>,
}
pub struct EcFftCosetPrecomputation<F: PrimeField, P: EcFftParameters<F>> {
    pub coset: Vec<F>,
    pub steps: Vec<EcFftPrecomputationStep<F, P>>,
}

pub struct EcFftPrecomputation<F: PrimeField, P: EcFftParameters<F>> {
    pub coset_precomputations: Vec<EcFftCosetPrecomputation<F, P>>,
}

impl<F: PrimeField, P: EcFftParameters<F>> EcFftCosetPrecomputation<F, P> {
    /// From `evals` the evaluations of a polynomial on `self.steps[0].s`,
    /// return the evaluations of the polynomial on `self.steps[0].s_prime` in `O(nlogn)`.
    /// See https://solvable.group/posts/ecfft/ for a simple explanation of this function.
    pub fn extend(&self, evals: &[F]) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }
        assert_eq!(
            n & (n - 1),
            0,
            "The number of evaluations should be a power of 2."
        );
        let log_n = n.trailing_zeros() as usize;
        assert!(
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
        } = &self.steps[self.steps.len() - log_n];
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
}

impl<F: PrimeField, P: EcFftParameters<F>> EcFftPrecomputation<F, P> {
    /// Evaluates polynomial of degree `<n` on the sub-coset of size `n` in O(nlog^2n).
    /// Expects the polynomial to have a power of two coefficients, so one may need to resize with zeros before calling this.
    pub fn evaluate_over_domain(&self, poly: &DensePolynomial<F>) -> Vec<F> {
        let n = poly.len();
        if n == 1 {
            return vec![poly.coeffs[0]];
        }
        assert_eq!(
            n & (n - 1),
            0,
            "The number of coefficients should be a power of 2."
        );
        let log_n = n.trailing_zeros() as usize;
        assert!(
            log_n <= P::LOG_N,
            "The polynomial can have degree at most {}.",
            1 << P::LOG_N
        );
        let precomputations = &self.coset_precomputations;
        let low = poly.coeffs[..n / 2].to_vec();
        let high = poly.coeffs[n / 2..].to_vec();
        let low_evals = self.evaluate_over_domain(&DensePolynomial { coeffs: low });
        let high_evals = self.evaluate_over_domain(&DensePolynomial { coeffs: high });
        let low_evals_prime = precomputations[P::LOG_N - log_n].extend(&low_evals);
        let high_evals_prime = precomputations[P::LOG_N - log_n].extend(&high_evals);

        let coset = &precomputations[P::LOG_N - log_n].coset;
        assert_eq!(n, coset.len());
        let mut ans = Vec::new();
        for i in 0..n / 2 {
            ans.push(low_evals[i] + coset[2 * i].pow([n as u64 / 2]) * high_evals[i]);
            ans.push(
                low_evals_prime[i] + coset[2 * i + 1].pow([n as u64 / 2]) * high_evals_prime[i],
            );
        }

        ans
    }
}
