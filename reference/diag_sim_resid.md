# Quantile Residuals of Simulated Responses

Draws `n_sim` responses with
[`simulate.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/simulate.StatmodFit.md)
and returns the quantile residual of each under the fitted parameters.

## Usage

``` r
diag_sim_resid(fit, n_sim, seed)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- n_sim:

  The number of replicates.

- seed:

  A single number or `NULL`.

## Value

A matrix with one row per observation and `n_sim` columns.
