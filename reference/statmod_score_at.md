# The Score of the Weighted Log-Likelihood

One block per distribution parameter, \\X_k'(w\\g_k)\\ with \\g_k\\ the
per-observation derivative of the log-density in the link-scale
predictor.

## Usage

``` r
statmod_score_at(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

## Value

A named list of gradient vectors, one per parameter.
