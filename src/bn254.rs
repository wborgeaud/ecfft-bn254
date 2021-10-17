use std::{convert::TryInto, marker::PhantomData};

use crate::{
    ecfft::{EcFftCosetPrecomputation, EcFftParameters, EcFftPrecomputationStep},
    utils::{isogeny::Isogeny, matrix::Matrix},
};
use ark_ff::{BigInteger256, Field};

pub type F = ark_bn254::Fq;
pub type FParams = ark_bn254::FqParameters;

pub struct Bn254EcFftParameters;

impl EcFftParameters<F> for Bn254EcFftParameters {
    const LOG_N: usize = 14;

    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F> {
        std::fs::read_to_string("bn254_coset")
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect::<Vec<u64>>()
            .chunks(4)
            .map(|chunk| BigInteger256::new(chunk.try_into().unwrap()).into())
            .collect()
    }

    fn isogenies() -> Vec<Isogeny<F>> {
        std::fs::read_to_string("bn254_isogenies")
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect::<Vec<u64>>()
            .chunks(5 * 4)
            .map(|chunk| {
                let numerator = (0..3)
                    .map(|i| {
                        BigInteger256::new(chunk[i * 4..(i + 1) * 4].try_into().unwrap()).into()
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                let denominator = (3..5)
                    .map(|i| {
                        BigInteger256::new(chunk[i * 4..(i + 1) * 4].try_into().unwrap()).into()
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                Isogeny {
                    numerator,
                    denominator,
                }
            })
            .collect()
    }

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
}

#[cfg(test)]
mod tests {
    use crate::ecfft::EcFftParameters;

    use crate::bn254::{Bn254EcFftParameters, F};
    use ark_ff::PrimeField;
    use ark_poly::{univariate::DensePolynomial, Polynomial};
    use ark_std::{
        rand::{distributions::Standard, prelude::Distribution, Rng},
        test_rng,
    };

    use super::*;

    #[test]
    fn test_precompute() {
        Bn254EcFftParameters::precompute_on_coset(&Bn254EcFftParameters::coset());
        Bn254EcFftParameters::precompute_on_coset(
            &Bn254EcFftParameters::coset()
                .into_iter()
                .step_by(2)
                .collect::<Vec<_>>(),
        );
    }

    fn test_extend_i<F: PrimeField, P: EcFftParameters<F>>(
        i: usize,
        precomputation: &EcFftCosetPrecomputation<F, P>,
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
        let precomputation =
            Bn254EcFftParameters::precompute_on_coset(&Bn254EcFftParameters::coset());
        for i in 1..Bn254EcFftParameters::LOG_N {
            test_extend_i::<F, _>(i, &precomputation);
        }
    }

    #[test]
    fn test_eval() {
        type P = Bn254EcFftParameters;
        let precomputation = P::precompute();
        for i in 0..P::LOG_N {
            let mut rng = test_rng();
            let coeffs: Vec<F> = (0..P::N >> i).map(|_| rng.gen()).collect();
            let poly = DensePolynomial { coeffs };
            let now = std::time::Instant::now();
            let evals = P::sub_coset(i)
                .iter()
                .map(|x| poly.evaluate(x))
                .collect::<Vec<_>>();
            dbg!(now.elapsed().as_secs_f32());
            assert_eq!(evals, precomputation.evaluate_over_domain(&poly));
            dbg!(now.elapsed().as_secs_f32());
        }
    }
}
