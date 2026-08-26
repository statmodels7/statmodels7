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
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- ep:

  The linear predictors and the parameters they imply, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- p:

  Which distribution parameter's equation, a string.

- expected:

  `TRUE` for the expected information, `FALSE` for the observed one.

- approx:

  How the expected information is approximated for a family with no
  closed form.

## Value

A list with `w` and `z`, each a numeric vector of length `spec@n_obs`.
`NULL` where the curvature is not usable, which is any non-finite or
non-positive \\h_i\\; the observed information can produce both far from
the optimum.

## Details

For a Gaussian response on the identity link the quadratic is exact and
one pass answers the problem. Elsewhere it is the local approximation a
scoring step works on, and the weights are rebuilt at each iteration.
