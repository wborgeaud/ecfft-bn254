use ark_ff::PrimeField;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial,
};

/// Returns the GCD of `a` and `b`.
pub fn euclids_algorithm<F: PrimeField>(
    a: DenseOrSparsePolynomial<F>,
    b: DenseOrSparsePolynomial<F>,
) -> DensePolynomial<F> {
    let mut r0 = a;
    let mut r1 = b;
    let mut r2 = DenseOrSparsePolynomial::from(r0.divide_with_q_and_r(&r1).unwrap().1);
    while !r2.is_zero() {
        r0 = r1;
        r1 = r2;
        r2 = DenseOrSparsePolynomial::from(r0.divide_with_q_and_r(&r1).unwrap().1);
    }
    r1.into()
}

/// Returns the extended GCD of `a` and `b`. That is polynomials `s` and `t` such that
/// `a * t + b * s = gcd(a, b)`.
pub fn extended_gcd_algorithm<F: PrimeField>(
    a: DenseOrSparsePolynomial<F>,
    b: DenseOrSparsePolynomial<F>,
) -> (DensePolynomial<F>, DensePolynomial<F>) {
    let mut s = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let mut old_s = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    let mut t = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    let mut old_t = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let mut r = b;
    let mut old_r = a;
    while !r.is_zero() {
        let (quotient, new_r) = old_r.divide_with_q_and_r(&r).unwrap();
        old_r = r;
        r = DenseOrSparsePolynomial::from(new_r);
        let new_s = &old_s - &(&quotient * &s);
        old_s = s;
        s = new_s;
        let new_t = &old_t - &(&quotient * &t);
        old_t = t;
        t = new_t;
    }
    (old_s, old_t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    use ark_ff::Zero;
    use ark_poly::univariate::DenseOrSparsePolynomial;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    #[test]
    fn test_euclids_algorithm() {
        let mut rng = test_rng();
        for _ in 0..100 {
            let a = DensePolynomial::<Fr>::rand(10, &mut rng);
            let b = DensePolynomial::<Fr>::rand(10, &mut rng);
            let gcd = euclids_algorithm::<Fr>(
                DenseOrSparsePolynomial::<Fr>::from(a.clone()),
                DenseOrSparsePolynomial::<Fr>::from(b.clone()),
            );
            let a_dors = DenseOrSparsePolynomial::<Fr>::from(a);
            let b_dors = DenseOrSparsePolynomial::<Fr>::from(b);
            let gcd_dors = DenseOrSparsePolynomial::<Fr>::from(gcd.clone());
            assert_eq!(
                a_dors.divide_with_q_and_r(&gcd_dors).unwrap().1,
                DensePolynomial::<Fr>::zero()
            );
            assert_eq!(
                b_dors.divide_with_q_and_r(&gcd_dors).unwrap().1,
                DensePolynomial::<Fr>::zero()
            );
        }
    }

    #[test]
    fn test_extended_gcd() {
        let mut rng = test_rng();
        for _ in 0..100 {
            let a = DensePolynomial::<Fr>::rand(10, &mut rng);
            let b = DensePolynomial::<Fr>::rand(10, &mut rng);
            let (s, t) = extended_gcd_algorithm::<Fr>(
                DenseOrSparsePolynomial::<Fr>::from(a.clone()),
                DenseOrSparsePolynomial::<Fr>::from(b.clone()),
            );
            let gcd = euclids_algorithm::<Fr>(
                DenseOrSparsePolynomial::<Fr>::from(a.clone()),
                DenseOrSparsePolynomial::<Fr>::from(b.clone()),
            );
            assert_eq!(&(&a * &t) + &(&b * &s), gcd);
        }
    }
}
