# ECFFT algorithms on the BN254 base field
This crate implements structs and traits to implement the ECFFT algorithms from the paper [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.

A concrete implementation is provided for the BN254 base field which is not FFT friendly (two-adicity of 1).

### References
- [Elliptic Curve Fast Fourier Transform (ECFFT) Part I: Fast Polynomial Algorithms over all Finite Fields](https://arxiv.org/abs/2107.08473) by Eli Ben-Sasson, Dan Carmon, Swastik Kopparty and David Levit.
- [The ECFFT algorithm](https://solvable.group/posts/ecfft/).
- [ECFFT on the BN254 base field in Rust](https://solvable.group/posts/ecfft-bn254/).
