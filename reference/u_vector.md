# The Trace of the Determinant's Movement With the Coefficients

\\u_c = \mathrm{tr}(M\\\partial K/\partial\beta_c)\\, assembled one
crossprod per distribution parameter.

## Usage

``` r
u_vector(
  spec,
  design,
  coef,
  M,
  params,
  npar,
  offs,
  total,
  d3 = NULL,
  G = NULL,
  key = NULL,
  rows = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- M:

  The matrix the trace is taken against.

- params:

  The distribution's parameter names.

- npar, offs, total:

  The block sizes, their offsets and the total.

- d3:

  The array \\\partial K/\partial\eta\\ is read from, or `NULL` for the
  family's third derivative at the fitted predictors.

- G:

  The per-observation diagonals of \\M\\, or `NULL`.

- key:

  The builder of `d3`'s component names, or `NULL`.

- rows:

  `NULL` for the contraction against each equation's design, or one
  matrix per distribution parameter giving the derivative of that
  equation's predictor in a wider vector, as
  [`filter_joint_movement()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_movement.md)
  returns them. The result is then as long as their common width.

## Value

A numeric vector as long as the stacked coefficients, or as the width of
`rows`.
