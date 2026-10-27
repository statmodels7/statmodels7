# The Penalized Information at a Point

Returns \\H + S\\ as a dense matrix and sets the penalty's non-finite
entries to zero.

## Usage

``` r
statmod_penalized_at(
  spec,
  coef,
  design,
  hyper,
  expected = TRUE,
  emethod = NULL
)
```

## Arguments

- spec:

  The specification.

- coef:

  A named list of coefficients, one vector per distribution parameter.

- design:

  The design.

- hyper:

  The hyperparameters.

- expected:

  Whether the expected information is used.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

A dense symmetric matrix over the coefficients of every equation.

## Details

This is the matrix from which a penalized fit's curvature is read: the
model's information from
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
and the penalty's Hessian from
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md).
It is written once because its readers, the rank test of
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
in
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
the readings of the alternation
([`alternation_readings()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md))
and the flat-direction test of the outer search, must use the same
matrix at the same point.

A kinked penalty contributes no curvature away from its kink, and any
non-finite entry would be the kink itself reached by a hair, so the
penalty passes through
[`zap_nonfinite()`](https://statmodels7.github.io/statmodels7/reference/zap_nonfinite.md)
as it does everywhere else.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md),
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
