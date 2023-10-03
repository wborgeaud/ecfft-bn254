use ark_bn254::Fr;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{rand::Rng, test_rng, time::Duration};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ecfft_bn254::bn254_scalar::{Bn254ScalarEcFftParameters, F};
use ecfft_bn254::ecfft::EcFftParameters;

fn interpolation(c: &mut Criterion) {
    type P = Bn254ScalarEcFftParameters;
    let precomputation = P::precompute();
    let mut rng = test_rng();

    let mut group = c.benchmark_group("interpolation");
    group.measurement_time(Duration::from_secs(30));

    for log_n in 1..=P::LOG_N {
        group.bench_with_input(BenchmarkId::new("ECFFT", log_n), &log_n, |b, _| {
            let evals: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            b.iter(|| precomputation.interpolate(&evals));
        });

        group.bench_with_input(BenchmarkId::new("Classic", log_n), &log_n, |b, _| {
            let evals: Vec<Fr> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let domain = Radix2EvaluationDomain::<F>::new(1 << log_n).unwrap();
            b.iter(|| domain.ifft(&evals));
        });
    }
}

criterion_group!(benches, interpolation);
criterion_main!(benches);
