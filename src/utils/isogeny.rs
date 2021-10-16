use ark_ff::PrimeField;

#[derive(Debug, Clone, Copy)]
pub struct Isogeny<F: PrimeField> {
    pub numerator: [F; 3],
    pub denominator: [F; 2],
}

impl<F: PrimeField> Isogeny<F> {
    pub fn eval(&self, x: F) -> F {
        let Isogeny {
            numerator: [a0, a1, a2],
            denominator: [b0, b1],
        } = self;
        (*a0 + *a1 * x + *a2 * x * x) / (*b0 + *b1 * x)
    }

    pub fn eval_den(&self, x: F) -> F {
        let Isogeny {
            denominator: [b0, b1],
            ..
        } = self;
        *b0 + *b1 * x
    }
}
