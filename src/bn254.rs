use std::{convert::TryInto, marker::PhantomData};

use crate::{
    ecfft::{EcFftParameters, EcFftPrecomputation, EcFftPrecomputationStep},
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

    fn precompute(coset: Vec<F>) -> EcFftPrecomputation<F, Self> {
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

        EcFftPrecomputation {
            steps,
            final_s: s[0],
            final_s_prime: s_prime[0],
            _phantom: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ecfft::EcFftParameters;

    use super::Bn254EcFftParameters;

    #[test]
    fn test_precompute() {
        Bn254EcFftParameters::precompute(Bn254EcFftParameters::coset());
        Bn254EcFftParameters::precompute(
            Bn254EcFftParameters::coset()
                .into_iter()
                .step_by(2)
                .collect(),
        );
    }
}
