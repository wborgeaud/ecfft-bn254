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
            let mut matrices_s_prime = Vec::new();
            let mut inverse_matrices_s = Vec::new();
            let mut matrices_s = Vec::new();
            let mut inverse_matrices_s_prime = Vec::new();

            for j in 0..nn {
                let (s0, s1) = (s[j], s[j + nn]);
                debug_assert_eq!(psi.eval(s0), psi.eval(s1));
                inverse_matrices_s.push(
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

                matrices_s.push(Matrix([
                    [
                        psi.eval_den(s0).pow([q as u64]),
                        s0 * psi.eval_den(s0).pow([q as u64]),
                    ],
                    [
                        psi.eval_den(s1).pow([q as u64]),
                        s1 * psi.eval_den(s1).pow([q as u64]),
                    ],
                ]));

                let (s0, s1) = (s_prime[j], s_prime[j + nn]);
                debug_assert_eq!(psi.eval(s0), psi.eval(s1));
                inverse_matrices_s_prime.push(
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

                matrices_s_prime.push(Matrix([
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
                matrices_s,
                inverse_matrices_s,
                matrices_s_prime,
                inverse_matrices_s_prime,
                vanish_on_s_prime: Vec::new(),
                c_evals: Vec::new(),
                _phantom: PhantomData,
            });
            s = s.into_iter().take(nn).map(|x| psi.eval(x)).collect();
            s_prime = s_prime.into_iter().take(nn).map(|x| psi.eval(x)).collect();
        }
        debug_assert_eq!((s.len(), s_prime.len()), (1, 1));

        let mut result = EcFftCosetPrecomputation {
            coset: coset.to_vec(),
            steps,
        };

        result.calc_vanish_on_s_prime();

        result
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
    pub matrices_s: Vec<Matrix<F>>,
    pub inverse_matrices_s: Vec<Matrix<F>>,
    pub matrices_s_prime: Vec<Matrix<F>>,
    pub inverse_matrices_s_prime: Vec<Matrix<F>>,
    pub vanish_on_s_prime: Vec<F>,
    pub c_evals: Vec<F>,
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
    /// return the evaluations of the polynomial on `self.steps[0].s_prime` in `O(n * log n)`.
    /// See https://solvable.group/posts/ecfft/ for a simple explanation of this function.
    pub fn extend_s_to_s_prime(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.extend_in_place_s_to_s_prime(&mut evals);
        evals
    }

    /// Mutate `evals`, which contains the evaluations of a polynomial on `self.steps[0].s`,
    /// to store the evaluations of the polynomial on `self.steps[0].s_prime` in `O(n * log n)`.
    /// See https://solvable.group/posts/ecfft/ for a simple explanation of this function.
    pub fn extend_in_place_s_to_s_prime(&self, evals: &mut [F]) {
        let n = evals.len();
        if n == 1 {
            return;
        }
        assert_eq!(
            n.next_power_of_two(),
            n,
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
            matrices_s_prime,
            inverse_matrices_s,
            ..
        } = &self.steps[self.steps.len() - log_n];
        let nn = n / 2;
        let (p0_evals, p1_evals) = evals.split_at_mut(nn);

        for (m, (p0, p1)) in inverse_matrices_s
            .iter()
            .zip(p0_evals.iter_mut().zip(p1_evals.iter_mut()))
        {
            m.multiply_in_place(p0, p1);
        }
        self.extend_in_place_s_to_s_prime(p0_evals);
        self.extend_in_place_s_to_s_prime(p1_evals);

        for (m, (p0, p1)) in matrices_s_prime
            .iter()
            .zip(p0_evals.iter_mut().zip(p1_evals.iter_mut()))
        {
            m.multiply_in_place(p0, p1)
        }
    }

    /// From `evals` the evaluations of a polynomial on `self.steps[0].s_prime`,
    /// return the evaluations of the polynomial on `self.steps[0].s` in `O(n * log n)`.
    /// See https://solvable.group/posts/ecfft/ for a simple explanation of this function.
    pub fn extend_s_prime_to_s(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.extend_in_place_s_prime_to_s(&mut evals);
        evals
    }

    /// Mutate `evals`, which contains the evaluations of a polynomial on `self.steps[0].s_prime`,
    /// to store the evaluations of the polynomial on `self.steps[0].s` in `O(n * log n)`.
    /// See https://solvable.group/posts/ecfft/ for a simple explanation of this function.
    pub fn extend_in_place_s_prime_to_s(&self, evals: &mut [F]) {
        let n = evals.len();
        if n == 1 {
            return;
        }
        assert_eq!(
            n.next_power_of_two(),
            n,
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
            matrices_s,
            inverse_matrices_s_prime,
            ..
        } = &self.steps[self.steps.len() - log_n];
        let nn = n / 2;
        let (p0_evals, p1_evals) = evals.split_at_mut(nn);

        for (m, (p0, p1)) in inverse_matrices_s_prime
            .iter()
            .zip(p0_evals.iter_mut().zip(p1_evals.iter_mut()))
        {
            m.multiply_in_place(p0, p1);
        }
        self.extend_in_place_s_prime_to_s(p0_evals);
        self.extend_in_place_s_prime_to_s(p1_evals);

        for (m, (p0, p1)) in matrices_s
            .iter()
            .zip(p0_evals.iter_mut().zip(p1_evals.iter_mut()))
        {
            m.multiply_in_place(p0, p1)
        }
    }

    fn calc_vanish_on_s_prime(&mut self) {
        for step in self.steps.iter_mut() {
            step.vanish_on_s_prime = step
                .s_prime
                .iter()
                .map(|val| step.s.iter().fold(F::one(), |acc, s| acc * (*val - *s)))
                .collect::<Vec<F>>();
        }
    }
}

impl<F: PrimeField, P: EcFftParameters<F>> EcFftPrecomputation<F, P> {
    /// Evaluates polynomial of degree `<n` on the sub-coset of size `n` in O(n * log^2 n).
    /// Expects the polynomial to have a power of two coefficients, so one may need to resize with zeros before calling this.
    pub fn evaluate_over_domain(&self, poly: &DensePolynomial<F>) -> Vec<F> {
        let mut evaluations = poly.to_vec();
        let mut scratch1 = poly.coeffs.clone();
        self.ecfft_in_place(&mut evaluations, &mut scratch1);
        evaluations
    }

    /// Evaluates polynomial of degree `<n` on the sub-coset of size `n` in O(n * log^2 n).
    /// Expects the polynomial to have a power of two coefficients, so one may need to resize with zeros before calling this.
    pub fn ecfft_in_place(&self, poly: &mut [F], scratch1: &mut [F]) {
        let n = poly.len();
        if n == 1 {
            return;
        }
        assert_eq!(
            n.next_power_of_two(),
            n,
            "The number of coefficients should be a power of 2."
        );

        let log_n = n.trailing_zeros() as usize;
        assert!(
            log_n <= P::LOG_N,
            "The polynomial can have degree at most {}.",
            1 << P::LOG_N
        );
        let precomputations = &self.coset_precomputations;
        let (low, high) = poly.split_at_mut(n / 2);
        let (low_1, high_1) = scratch1.split_at_mut(n / 2);
        self.ecfft_in_place(low, low_1);
        self.ecfft_in_place(high, high_1);
        low_1.copy_from_slice(low);
        high_1.copy_from_slice(high);
        let coset = &precomputations[P::LOG_N - log_n].coset;
        assert_eq!(n, coset.len());
        (0..(n / 2)).for_each(|i| {
            poly[2 * i] = low_1[i] + coset[2 * i].pow([n as u64 / 2]) * high_1[i];
        });

        precomputations[P::LOG_N - log_n].extend_in_place_s_to_s_prime(low_1);
        precomputations[P::LOG_N - log_n].extend_in_place_s_to_s_prime(high_1);

        (0..(n / 2)).for_each(|i| {
            poly[2 * i + 1] = low_1[i] + coset[2 * i + 1].pow([n as u64 / 2]) * high_1[i];
        });
    }
}

impl<F: PrimeField, P: EcFftParameters<F>> EcFftPrecomputation<F, P> {
    /// Performs the REDC algorithm from https://arxiv.org/pdf/2107.08473.pdf
    /// Computes the evaluations of Q(x) on coset, where Q(x) is congruent to
    /// P(x)*Z_0(x)^-1 mod X^n/2. We assume evals are input in the same manner as s and s_prime
    /// are calculated.
    pub fn redc_in_place(&self, evals: &mut [F]) {
        let n = evals.len();

        // If there is only one evaluation then the polynomial is constant and so reduction by X^n/2 does nothing.
        if n == 1 {
            return;
        }

        assert_eq!(
            n.next_power_of_two(),
            n,
            "The number of coefficients should be a power of 2."
        );

        let log_n = n.trailing_zeros() as usize;
        assert!(
            log_n <= P::LOG_N,
            "The polynomial can have degree at most {}.",
            1 << P::LOG_N
        );
        let precomputations = &self.coset_precomputations;

        let coset = &precomputations[P::LOG_N - log_n].coset;
        let vanish_on_s_prime = if let Some(step) = precomputations[P::LOG_N - log_n].steps.first()
        {
            step.vanish_on_s_prime.clone()
        } else {
            vec![coset[1] - coset[0]]
        };
        assert_eq!(n, coset.len());

        let mut evals_s = evals.iter().step_by(2).copied().collect::<Vec<F>>();

        let mut evals_s_prime = evals.iter().skip(1).step_by(2).copied().collect::<Vec<F>>();

        (0..(n / 2)).for_each(|i| {
            evals_s[i] = coset[2 * i].pow([n as u64 / 2]).inverse().unwrap() * evals_s[i];
        });

        let g_on_s_prime = precomputations[P::LOG_N - log_n].extend_s_to_s_prime(&evals_s);

        (0..(n / 2)).for_each(|i| {
            evals_s_prime[i] = (evals_s_prime[i]
                - g_on_s_prime[i] * coset[2 * i + 1].pow([n as u64 / 2]))
                * vanish_on_s_prime[i].inverse().unwrap();
        });

        let h_0 = &precomputations[P::LOG_N - log_n].extend_s_prime_to_s(&evals_s_prime);

        for (i, eval) in evals.iter_mut().enumerate() {
            if i % 2 == 0 {
                *eval = h_0[i / 2];
            } else {
                *eval = evals_s_prime[i / 2];
            }
        }
    }

    pub fn redc(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.redc_in_place(&mut evals);
        evals
    }
}
