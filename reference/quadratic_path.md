# The Path of a Quadratic Penalty in Closed Form

The Path of a Quadratic Penalty in Closed Form

## Usage

``` r
quadratic_path(fit, sys, u, lam_hat, log_lambda)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- sys:

  The list from
  [`penalty_system()`](https://statmodels7.github.io/statmodels7/reference/penalty_system.md).

- u:

  The penalized unit of the term.

- lam_hat:

  The fitted penalty parameter.

- log_lambda:

  `NULL` or the grid of \\\log\lambda\\.

## Value

A list with `log_lambda` and `coef`, as
[`kinked_path()`](https://statmodels7.github.io/statmodels7/reference/kinked_path.md).
