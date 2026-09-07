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
It is written once because three readers want the same matrix at the
same point and must not disagree about it –
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)'s
smoother,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)'s
variance and
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)'s
rank test.

A kinked penalty contributes no curvature away from its kink, and any
non-finite entry would be the kink itself reached by a hair, so the
penalty passes through
[`zap_nonfinite()`](https://statmodels7.github.io/statmodels7/reference/zap_nonfinite.md)
as it does everywhere else.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md),
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
