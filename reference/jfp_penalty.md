# The Jeffreys Penalty of a Model, With the Sign of a Penalty

Computes \\-\sum_b a_b\log\|\tilde I_b\|\\, the quantity that the
objective adds to the negative log-likelihood, or its gradient or
Hessian over the stacked coefficients.

## Usage

``` r
jfp_penalty(
  spec,
  design,
  coef,
  what = c("value", "gradient", "hessian", "parts"),
  keys = NULL,
  joint = FALSE
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- what:

  `"value"`, `"gradient"`, `"hessian"` or `"parts"`.

- keys:

  As in
  [`jfp_blocks()`](https://statmodels7.github.io/statmodels7/reference/jfp_blocks.md).

- joint:

  Whether to return the whole joint vector beside a filter.

## Value

`0`, a zero vector, `NULL` or `NULL` where the model carries no Jeffreys
term; otherwise a number, a numeric vector over the stacked
coefficients, a base matrix, or the list of
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md).

## Details

Where a block's information is not positive definite the value is `Inf`,
so a line search rejects the point, and the gradient and the Hessian
signal an error.

`"parts"` returns the Hessian split as
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
gives it, with the sign of the penalty: the objective's curvature gains
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
at `H = -weights`, plus `-refresh`, plus `tcrossprod(lowrank)`.

Beside a filter the state spans the joint vector of the coefficients and
the filter's free parameters. With `joint = FALSE` the result is then
restricted to the coefficients, which is what the objective over the
coefficients reads;
[`jfp_joint_extra()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_extra.md)
returns the rest.
