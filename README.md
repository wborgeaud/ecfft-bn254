# ECFFT algorithms on the BN254 base field

This crate implements structs and traits for the ECFFT algorithms from the paper [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.

A concrete implementation is provided for the BN254 base field which is not FFT friendly (two-adicity of 1).

### Example

```rust
fn test_evaluations() {
    type P = Bn254EcFftParameters;
    // ECFFT precomputations.
    let precomputation = P::precompute();
    // Can interpolate polynomials up to degree 2^14.
    let log_n = 14;
    let mut rng = test_rng();
    // Generate a random polynomial.
    let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
    let poly = DensePolynomial { coeffs };
    // Naive evaluations.
    let evals = P::coset()
        .iter()
        .map(|x| poly.evaluate(x))
        .collect::<Vec<_>>();
    // ECFFT evaluations.
    let ecfft_evals = precomputation.evaluate_over_domain(&poly);

    assert_eq!(evals, ecfft_evals);
}
```

### BLS12-381

The base field of the BLS12-381 curve is also supported, also for degrees up to `2^14`. Credits to [Youssef El Housni](https://github.com/yelhousni) for [finding a curve with 2-adicity 14](https://ethresear.ch/t/bw6-over-bls12-381/10321/4).

### Precomputations

The implementation uses precomputations for the coset and isogenies used in the ECFFT. These precomputations are computed in `get_params.sage` and are stored in the `bn254_coset` and `bn254_isogenies` files.

To implement the ECFFT for other fields, similar precomputations should be performed. For example, here is how to generate the precomputations for BLS12-381:

```bash
# sage get_params.sage p a b output_filename
sage get_params.sage 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab 0x1800fb41dab7368489a980e14a746abfe7c87588aac25c113301d524b734a5043bbc89dd7d0c5b41de5d348ac2e838c6 0x11c65a0a6e52b8b88366e0b0df28c6804f14f35cb833cb0d918c9e758f044d95777beb965a967af4ef518ad0618a809a bls12-381
```

### Benchmarks

Here is a comparison of the running time for the evaluation of a polynomial of degree `n-1` on a domain of `n` points using 3 algorithms:

- the naive evaluation in `O(n^2)`,
- the classic FFT (on the FFT-friendly BN254 scalar field) in `O(n * log n)`,
- the ECFFT ENTER algorithm in `O(n * log^2 n)`.

| `log n` | Naive (ms)  | Classic (ms) | ECFFT (ms) | Naive/ECFFT | ECFFT/Classic |
| ------- | ----------- | ------------ | ---------- | ----------- | ------------- |
| 1       | 0.000165    | 0.000126     | 0.000384   | 0.429       | 3.056         |
| 2       | 0.00046     | 0.000256     | 0.002144   | 0.214       | 8.36          |
| 3       | 0.00203     | 0.000639     | 0.008599   | 0.236       | 13.456        |
| 4       | 0.00688     | 0.001781     | 0.030458   | 0.226       | 17.103        |
| 5       | 0.032354    | 0.003268     | 0.085556   | 0.378       | 26.177        |
| 6       | 0.119391    | 0.007594     | 0.239939   | 0.498       | 31.595        |
| 7       | 0.479542    | 0.018378     | 0.613242   | 0.782       | 33.368        |
| 8       | 1.873195    | 0.043694     | 1.425794   | 1.314       | 32.632        |
| 9       | 7.619662    | 0.093        | 3.964933   | 1.922       | 42.634        |
| 10      | 30.034845   | 0.20955      | 9.308925   | 3.226       | 44.423        |
| 11      | 121.564343  | 0.453727     | 22.186604  | 5.479       | 48.899        |
| 12      | 482.728362  | 0.976134     | 51.505625  | 9.372       | 52.765        |
| 13      | 1930.495799 | 2.166843     | 119.317395 | 16.18       | 55.065        |
| 14      | 7745.103265 | 4.57555      | 275.499648 | 28.113      | 60.211        |

### References

- [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.
- [The ECFFT algorithm](https://solvable.group/posts/ecfft/).
- [ECFFT on the BN254 base field in Rust](https://solvable.group/posts/ecfft-bn254/).
