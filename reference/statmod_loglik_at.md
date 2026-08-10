# The Weighted Log-Likelihood of a Specification at Given Coefficients

\\\sum_i w_i \log f(y_i; \theta_i)\\, the weights entering as given.

## Usage

``` r
statmod_loglik_at(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design; recomputed when absent.

## Value

A single number.
