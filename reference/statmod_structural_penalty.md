# The Penalty Over a Structural Term's Own Parameters

The value, gradient or Hessian of the penalties that a structural term
declares, in the term's own parameters.
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
evaluates the penalties on the coefficients.

## Usage

``` r
statmod_structural_penalty(
  spec,
  design,
  hyper,
  what = c("value", "gradient", "hessian")
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose structural state holds each term's current
  parameters.

- hyper:

  The hyperparameters, per penalized term.

- what:

  Which quantity: `"value"`, `"gradient"` or `"hessian"`.

## Value

A named list, one entry per structural term, each a numeric vector or
matrix over that term's parameters in their own order; empty when no
structural term carries a penalty.

## Details

The objective of the structural block includes these penalties through
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md),
so that its gradient has to include their derivative. Without it the
objective and its gradient describe different functions, and the
gradient check of optimizers7 reports that the objective changes at a
rate that the gradient does not predict.

The penalty is read on the unconstrained scale, which is the scale on
which the term carries its parameters.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
