# The Working Response and Weights of One Equation

\\h_i\\ and \\z_i = \eta_i + s_i/h_i\\, the weighted least squares
problem the log-likelihood is locally.

## Usage

``` r
coord_working(spec, ep, coef, design, p, expected, approx)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- ep:

  The linear predictors and parameters, from
  [`statmod_eta`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).

- coef:

  The coefficients.

- design:

  The design.

- p:

  Which distribution parameter.

- expected:

  Whether the information is the expected one.

- approx:

  How it is approximated.

## Value

A list with `w` and `z`, or `NULL` where the curvature is not usable.
