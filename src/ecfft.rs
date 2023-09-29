use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use num_bigint::BigUint;

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

            let vanish_on_s_prime = s_prime
                .iter()
                .map(|val| s.iter().fold(F::one(), |acc, s| acc * (*val - *s)))
                .collect::<Vec<F>>();

            let c_evals_s_prime = vanish_on_s_prime
                .iter()
                .enumerate()
                .map(|(i, val)| {
                    let field_val = val.pow([2 as u64]);
                    let big_int: BigUint = field_val.into();
                    let reduce_by: BigUint = (s[i].pow([nn as u64])).into();
                    let modulo = big_int % reduce_by;
                    F::from(modulo)
                })
                .collect::<Vec<F>>();

            let mut c_evals = Vec::new();
            let c_evals_s = vec![F::zero(); s.len()];

            for (i, c_eval) in c_evals.iter_mut().enumerate() {
                if i % 2 == 0 {
                    *c_eval = c_evals_s[i / 2];
                } else {
                    *c_eval = c_evals_s_prime[i / 2];
                }
            }
            steps.push(EcFftPrecomputationStep::<F, Self> {
                s: s.clone(),
                s_prime: s_prime.clone(),
                matrices_s,
                inverse_matrices_s,
                matrices_s_prime,
                inverse_matrices_s_prime,
                vanish_on_s_prime,
                c_evals,
                _phantom: PhantomData,
            });
            s = s.into_iter().take(nn).map(|x| psi.eval(x)).collect();
            s_prime = s_prime.into_iter().take(nn).map(|x| psi.eval(x)).collect();
        }
        debug_assert_eq!((s.len(), s_prime.len()), (1, 1));

        EcFftCosetPrecomputation {
            coset: coset.to_vec(),
            vanish_on_s_prime: Vec::default(),
            z0z0_rem_xnn: Vec::default(),
            z1z1_rem_xnn: Vec::default(),
            steps,
        }
    }

    /// Computes the ECFFT precomputations of all `Self::sub_coset(i)`.
    fn precompute() -> EcFftPrecomputation<F, Self> {
        let mut coset_precomputations = Vec::new();
        let coset = Self::coset();
        for i in 0..Self::LOG_N {
            let n = 2usize.pow(i as u32 + 1);

            let current_coset = coset
                .iter()
                .step_by(2usize.pow(Self::LOG_N as u32 - 1 - i as u32))
                .copied()
                .collect::<Vec<F>>();
            let mut precompute = Self::precompute_on_coset(&current_coset);
            if i == 0 {
                precompute.vanish_on_s_prime = vec![current_coset[1] - current_coset[0]];
                precompute.z0z0_rem_xnn = vec![current_coset[0].square(); 2];
                precompute.z1z1_rem_xnn = vec![current_coset[1].square(); 2];
            } else {
                let previous_precomp: &EcFftCosetPrecomputation<F, Self> =
                    &coset_precomputations[0];

                let s = current_coset.iter().step_by(2).copied().collect::<Vec<_>>();
                let s_prime = current_coset
                    .iter()
                    .skip(1)
                    .step_by(2)
                    .copied()
                    .collect::<Vec<_>>();

                precompute.vanish_on_s_prime = s_prime
                    .iter()
                    .map(|val| s.iter().fold(F::one(), |acc, s| acc * (*val - *s)))
                    .collect::<Vec<F>>();

                let z0_rem_xnnnn_sq_s0 = previous_precomp
                    .z0z0_rem_xnn
                    .iter()
                    .zip(previous_precomp.z1z1_rem_xnn.iter())
                    .map(|(s0, s1)| *s0 * *s1)
                    .collect::<Vec<F>>();
                let z0z0_rem_xnnnn_s0 = previous_precomp.modulo_xnn(z0_rem_xnnnn_sq_s0.as_slice());
                let z0z0_rem_xnnnn_s1 =
                    precompute.extend_s_to_s_prime(z0z0_rem_xnnnn_s0.as_slice());
                let z0z0_rem_xnnnn_coset = z0z0_rem_xnnnn_s0
                    .iter()
                    .zip(z0z0_rem_xnnnn_s1.iter())
                    .flat_map(|(s0, s1)| vec![*s0, *s1])
                    .collect::<Vec<F>>();
                let z0_coset = precompute
                    .vanish_on_s_prime
                    .iter()
                    .flat_map(|val| vec![F::zero(), *val])
                    .collect::<Vec<F>>();
                let z0_rem_xnn_coset_squared = z0_coset
                    .iter()
                    .zip(current_coset.iter())
                    .map(|(z0, x)| (*z0 - x.pow([n as u64 / 2])).square())
                    .collect::<Vec<F>>();
                let complicated_term = z0_rem_xnn_coset_squared
                    .iter()
                    .zip(z0z0_rem_xnnnn_coset.iter().zip(current_coset.iter()))
                    .map(|(z0_rem_xnn_squared, (z0z0_rem_xnnnn, x))| {
                        (*z0_rem_xnn_squared - *z0z0_rem_xnnnn)
                            * (x.pow([n as u64 / 4])).inverse().unwrap()
                    })
                    .collect::<Vec<F>>();
                let complex_term_two =
                    precompute.modulo_xnnnn(complicated_term.as_slice(), &z0z0_rem_xnnnn_coset);
                precompute.z0z0_rem_xnn = z0z0_rem_xnnnn_coset
                    .iter()
                    .zip(complex_term_two.iter().zip(current_coset.iter()))
                    .map(|(z0z0_rem_xnnnn, (complex_term_two, x))| {
                        *z0z0_rem_xnnnn + (x.pow([n as u64 / 4]) * complex_term_two)
                    })
                    .collect::<Vec<F>>();

                let vanish_s_prime_on_s = s
                    .iter()
                    .map(|val| s_prime.iter().fold(F::one(), |acc, s| acc * (*val - *s)))
                    .collect::<Vec<F>>();

                let z1_coset = vanish_s_prime_on_s
                    .iter()
                    .flat_map(|val| vec![*val, F::zero()])
                    .collect::<Vec<F>>();
                let z1z1 = z1_coset
                    .iter()
                    .zip(current_coset.iter())
                    .map(|(z1, x)| (*z1 - x.pow([n as u64 / 2])).square())
                    .collect::<Vec<F>>();
                precompute.z1z1_rem_xnn = precompute.modulo_xnn(&z1z1);
            }
            coset_precomputations.insert(0, precompute);
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
    pub vanish_on_s_prime: Vec<F>,
    pub z0z0_rem_xnn: Vec<F>,
    pub z1z1_rem_xnn: Vec<F>,
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

    /// Performs the REDC algorithm from https://arxiv.org/pdf/2107.08473.pdf
    /// Computes the evaluations of Q(x) on coset, where Q(x) is congruent to
    /// P(x)*Z_0(x)^-1 mod X^n/2. We assume evals are input in the same manner as s and s_prime
    /// are calculated.
    fn redc_in_place_xnn(&self, evals: &mut [F]) {
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

        let coset = &self.coset;
        let vanish_on_s_prime = &self.vanish_on_s_prime;
        assert_eq!(n, coset.len());

        let mut evals_s = evals.iter().step_by(2).copied().collect::<Vec<F>>();

        let mut evals_s_prime = evals.iter().skip(1).step_by(2).copied().collect::<Vec<F>>();

        (0..(n / 2)).for_each(|i| {
            evals_s[i] = coset[2 * i].pow([n as u64 / 2]).inverse().unwrap() * evals_s[i];
        });

        let g_on_s_prime = self.extend_s_to_s_prime(&evals_s);

        (0..(n / 2)).for_each(|i| {
            evals_s_prime[i] = (evals_s_prime[i]
                - g_on_s_prime[i] * coset[2 * i + 1].pow([n as u64 / 2]))
                * vanish_on_s_prime[i].inverse().unwrap();
        });

        let h_0 = self.extend_s_prime_to_s(&evals_s_prime);

        for (i, eval) in evals.iter_mut().enumerate() {
            if i % 2 == 0 {
                *eval = h_0[i / 2];
            } else {
                *eval = evals_s_prime[i / 2];
            }
        }
    }

    /// Performs the REDC algorithm from https://arxiv.org/pdf/2107.08473.pdf
    /// Computes the evaluations of Q(x) on coset, where Q(x) is congruent to
    /// P(x)*Z_0(x)^-1 mod X^n/2. We assume evals are input in the same manner as s and s_prime
    /// are calculated.
    fn redc_in_place_xnnnn(&self, evals: &mut [F]) {
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

        let coset = &self.coset;
        let vanish_on_s_prime = &self.vanish_on_s_prime;
        assert_eq!(n, coset.len());

        let mut evals_s = evals.iter().step_by(2).copied().collect::<Vec<F>>();

        let mut evals_s_prime = evals.iter().skip(1).step_by(2).copied().collect::<Vec<F>>();

        (0..(n / 2)).for_each(|i| {
            evals_s[i] = coset[2 * i].pow([n as u64 / 4]).inverse().unwrap() * evals_s[i];
        });

        let g_on_s_prime = self.extend_s_to_s_prime(&evals_s);

        (0..(n / 2)).for_each(|i| {
            evals_s_prime[i] = (evals_s_prime[i]
                - g_on_s_prime[i] * coset[2 * i + 1].pow([n as u64 / 4]))
                * vanish_on_s_prime[i].inverse().unwrap();
        });

        let h_0 = self.extend_s_prime_to_s(&evals_s_prime);

        for (i, eval) in evals.iter_mut().enumerate() {
            if i % 2 == 0 {
                *eval = h_0[i / 2];
            } else {
                *eval = evals_s_prime[i / 2];
            }
        }
    }

    fn modulo_xnn(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.modulo_in_place_xnn(&mut evals);
        evals
    }

    fn modulo_in_place_xnn(&self, evals: &mut [F]) {
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

        let c_evals = &self.z0z0_rem_xnn;

        self.redc_in_place_xnn(evals);

        evals
            .iter_mut()
            .zip(c_evals.iter())
            .for_each(|(val, c_eval)| *val = *val * *c_eval);

        self.redc_in_place_xnn(evals);
    }

    pub fn modulo_xnnnn(&self, evals: &[F], c_poly_evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.modulo_in_place_xnnnn(&mut evals, c_poly_evals);
        evals
    }

    fn modulo_in_place_xnnnn(&self, evals: &mut [F], c_poly_evals: &[F]) {
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

        self.redc_in_place_xnnnn(evals);

        evals
            .iter_mut()
            .zip(c_poly_evals.iter())
            .for_each(|(val, c_eval)| *val = *val * *c_eval);

        self.redc_in_place_xnnnn(evals);
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

    pub fn interpolate(&self, evals: &[F]) -> DensePolynomial<F> {
        let mut evals = evals.to_vec();
        self.interpolate_in_place(&mut evals);
        DensePolynomial::<F>::from_coefficients_vec(evals)
    }

    pub fn interpolate_in_place(&self, evals: &mut Vec<F>) {
        let n = evals.len();
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
        let s = coset.iter().step_by(2).copied().collect::<Vec<F>>();

        let u_on_coset = self.modulo(evals);

        let u_on_s = u_on_coset.iter().step_by(2).copied().collect::<Vec<F>>();
        let u_poly = self.interpolate(&u_on_s);
        let evals_on_s = evals.iter().step_by(2).copied().collect::<Vec<F>>();
        let v_on_s = evals_on_s
            .iter()
            .zip(u_on_s.iter().zip(s.iter()))
            .map(|(eval_on_s, (u_on_s, x))| {
                (*eval_on_s - *u_on_s) * (x.pow([n as u64 / 2])).inverse().unwrap()
            })
            .collect::<Vec<F>>();
        let v_poly = self.interpolate(&v_on_s);

        *evals = [u_poly.coeffs, v_poly.coeffs].concat();
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

        precomputations[P::LOG_N - log_n].redc_in_place_xnn(evals);
    }

    pub fn redc(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.redc_in_place(&mut evals);
        evals
    }

    pub fn modulo(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.modulo_in_place(&mut evals);
        evals
    }

    pub fn modulo_in_place(&self, evals: &mut [F]) {
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
        precomputations[P::LOG_N - log_n].modulo_in_place_xnn(evals);
    }
}
