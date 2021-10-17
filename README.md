# ECFFT algorithms on the BN254 base field
This crate implements structs and traits to implement the ECFFT algorithms from the paper [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.

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

### Precomputations
The implementation uses precomputations for the coset and isogenies used in the ECFFT. These precomputations are computed in `get_params.sage` and are stored in the `bn254_coset` and `bn254_isogenies` files.

To implement the ECFFT for other fields, similar precomputations should be performed.

### References
- [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.
- [The ECFFT algorithm](https://solvable.group/posts/ecfft/).
- [ECFFT on the BN254 base field in Rust](https://solvable.group/posts/ecfft-bn254/).
