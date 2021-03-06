use std::convert::TryInto;

use crate::{ecfft::EcFftParameters, utils::isogeny::Isogeny};
use ark_ff::BigInteger256;

pub type F = ark_bn254::Fq;
/// Number of 64-bit limbs needed to represent field elements.
const NUM_LIMBS: usize = 4;

/// ECFFT parameters for the BN254 base field `F`.
/// Computed with the curve `E = EllipticCurve(F, [a, b])` with
/// `a, b = 1, 5612291247948481584627780310922020304781354847659642188369727566000581075360`.
pub struct Bn254EcFftParameters;

impl EcFftParameters<F> for Bn254EcFftParameters {
    /// The curve `E` has order `21888242871839275222246405745257275088712935808829559400805562964428910444544`
    /// with factorization `2^14 * 3^2 * 229 * 503 * 205460939795467 * 55374745393148401254803 * 113267149255983544517087125127`.
    const LOG_N: usize = 14;

    const N: usize = 1 << Self::LOG_N;

    /// Get the coset from the `bn254_coset` file. This file can be generated by running `get_params.sage`.
    fn coset() -> Vec<F> {
        std::fs::read_to_string("bn254_coset")
            .expect("Run `get_params.sage` to generate the coset.")
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect::<Vec<u64>>()
            .chunks(NUM_LIMBS)
            .map(|chunk| BigInteger256::new(chunk.try_into().unwrap()).into())
            .collect()
    }

    /// Get the isogenies from the `bn254_isogenies` file. This file can be generated by running `get_params.sage`.
    fn isogenies() -> Vec<Isogeny<F>> {
        std::fs::read_to_string("bn254_isogenies")
            .expect("Run `get_params.sage` to generate the coset.")
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect::<Vec<u64>>()
            .chunks(5 * NUM_LIMBS)
            .map(|chunk| {
                let numerator = (0..3)
                    .map(|i| {
                        BigInteger256::new(
                            chunk[i * NUM_LIMBS..(i + 1) * NUM_LIMBS]
                                .try_into()
                                .unwrap(),
                        )
                        .into()
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                let denominator = (3..5)
                    .map(|i| {
                        BigInteger256::new(
                            chunk[i * NUM_LIMBS..(i + 1) * NUM_LIMBS]
                                .try_into()
                                .unwrap(),
                        )
                        .into()
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
}

#[cfg(test)]
mod tests {
    use crate::ecfft::{EcFftCosetPrecomputation, EcFftParameters, EcFftPrecomputationStep};

    use crate::bn254::{Bn254EcFftParameters, F};
    use ark_ff::PrimeField;
    use ark_poly::{univariate::DensePolynomial, Polynomial};
    use ark_std::{
        rand::{distributions::Standard, prelude::Distribution, Rng},
        test_rng,
    };

    #[test]
    /// Tests that precomputations don't panic.
    fn test_precompute() {
        Bn254EcFftParameters::precompute_on_coset(&Bn254EcFftParameters::coset());
        Bn254EcFftParameters::precompute_on_coset(
            &Bn254EcFftParameters::coset()
                .into_iter()
                .step_by(2)
                .collect::<Vec<_>>(),
        );
    }

    /// Tests the extend function with a polynomial of degree `2^i - 1`.
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
    /// Tests the extend function for various degrees.
    fn test_extend() {
        let precomputation =
            Bn254EcFftParameters::precompute_on_coset(&Bn254EcFftParameters::coset());
        for i in 1..Bn254EcFftParameters::LOG_N {
            test_extend_i::<F, _>(i, &precomputation);
        }
    }

    #[test]
    /// Tests the `evaluate_over_domain` function for various degrees.
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
