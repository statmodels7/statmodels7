# The Penalized Information at a Point

\\H + S\\, dense, with the penalty's non-finite entries zeroed.

## Usage

``` r
statmod_penalized_at(
  spec,
  coef,
  design,
  hyper,
  expected = TRUE,
  approx = "opg"
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

- approx:

  How the expected information is approximated.

## Value

A dense symmetric matrix over the coefficients of every equation.

## Details

The matrix a penalized fit's curvature is read from: the model's
information from
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
