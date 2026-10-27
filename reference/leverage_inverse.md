# The Inverse Penalized Information of the Hat Matrix

Computes \\K^{-1}\\, with \\K\\ the information of the observations plus
the Hessian of the penalty, over the coordinates `cols`. It is the
matrix that the hat matrix projects with, and not the Bayesian variance
of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which for a coordinate that the marginal criterion estimates holds the
curvature of that criterion.

## Usage

``` r
leverage_inverse(fit, design, cols, expected, emethod, labels = NULL)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- cols:

  The positions, among the stacked coefficients, of the coordinates that
  the variance covers.

- expected, emethod:

  As in
  [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md).

- labels:

  The names of the coordinates, for error messages.

## Value

A dense matrix.
