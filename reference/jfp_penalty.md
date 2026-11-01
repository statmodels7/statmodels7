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

The information \\\tilde I_b\\ is the family's expected information
without modification, as in the definition of Jeffreys. It generally
depends on the other parameters of the distribution as well, so the
value, the gradient and the Hessian have entries in the coefficients of
those parameters, and the outer criterion reads the value and its
derivatives in the coordinates that it estimates
([`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)).
The prior therefore also moves those coordinates, by an amount of order
\\1/n\\: for a gaussian with the block in the mean the REML estimate of
\\\sigma^2\\ is the residual sum of squares divided by \\n - p +
2a\tilde p\\, and in simulations over seven families the median shift of
a dispersion or shape parameter was 0.1 to 0.3 standard errors at \\n =
50\\ and 0.02 to 0.06 at \\n = 1000\\. A block in a parameter other than
the mean moves the mean only when the information of that parameter
depends on the mean.

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
