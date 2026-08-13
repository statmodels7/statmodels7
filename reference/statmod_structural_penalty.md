# The Penalty Over a Structural Term's Own Parameters

The derivative of the penalties a structural term declares, in the
term's own parameters rather than in the coefficients.

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
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- hyper:

  The hyperparameters.

- what:

  One of `"value"`, `"gradient"` or `"hessian"`.

## Value

A named list, one entry per structural term, each a numeric vector or
matrix over that term's parameters in their own order; empty when no
structural term carries a penalty.

## Details

The objective of the structural block includes these penalties through
[`statmod_penalty_at`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md),
so its gradient must include their derivative: without it the two
describe different functions, an optimizer walks until its budget runs
out, and `optimizers7`'s own check reports that the objective changes at
a rate the gradient does not predict.

The penalty is read on the UNCONSTRAINED scale, which is where the term
carries its parameters and where a deviation from a population value is
defined; for a deviation, whose link is the identity, the two scales
coincide.

## See also

[`statmod_penalty_at`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
