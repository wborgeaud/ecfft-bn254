use ark_ff::Field;

#[derive(Debug, Clone, Copy)]
/// X-coordinate of an isogeny of degree 2.
pub struct Isogeny<F: Field> {
    /// Coefficients of the isogeny's numerator.
    pub numerator: [F; 3],
    /// Coefficients of the isogeny's denominator.
    pub denominator: [F; 2],
}

impl<F: Field> Isogeny<F> {
    /// Evaluate the (x-coordinate of the) isogeny on a field element.
    pub fn eval(&self, x: F) -> F {
        let Isogeny {
            numerator: [a0, a1, a2],
            denominator: [b0, b1],
        } = self;
        (*a0 + *a1 * x + *a2 * x * x) / (*b0 + *b1 * x)
    }

    /// Evaluate isogeny's denominator on a field element.
    pub fn eval_den(&self, x: F) -> F {
        let Isogeny {
            denominator: [b0, b1],
            ..
        } = self;
        *b0 + *b1 * x
    }
}
