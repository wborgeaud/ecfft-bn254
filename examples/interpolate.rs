use ark_std::{rand::Rng, test_rng};
use ecfft_bn254::bn254_scalar::{Bn254ScalarEcFftParameters, F};
use ecfft_bn254::ecfft::EcFftParameters;

fn main() {
    type P = Bn254ScalarEcFftParameters;
    let precomputation = P::precompute();
    let mut rng = test_rng();

    let log_n = P::LOG_N;
    let evals: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
    precomputation.interpolate(&evals);
}
