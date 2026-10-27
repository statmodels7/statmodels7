# The Path of a Kinked Penalty by Refitting

Fits the model at each value of \\\lambda\\ with the other
hyperparameters and the coefficients of the other equations held at
their fitted values, from the largest value to the smallest and each fit
from the previous one.

## Usage

``` r
kinked_path(fit, sys, u, lam_hat, log_lambda)
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

A list with `log_lambda` (increasing) and `coef`, a matrix with one row
per value and one column per drawn coefficient.
