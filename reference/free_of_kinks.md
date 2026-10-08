# A Score Over the Coordinates a Kink Leaves Free

Zeroes the entries of a score at the coordinates
[`zero_kinked()`](https://statmodels7.github.io/statmodels7/reference/zero_kinked.md)
returns.

## Usage

``` r
free_of_kinks(score, spec, design, coef)
```

## Arguments

- score:

  A score over the stacked coefficients.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- coef:

  The coefficients, a named list.

## Value

`score`, with those entries set to zero.

## Details

At a coordinate a kinked penalty holds at zero the objective has no
derivative, and the entry the smooth part reports is the
log-likelihood's score, which the kink's subgradient interval contains
rather than cancels. Read as a residual it inflated the mode error of a
fit sitting at its mode: on a lasso over a dispersion's ten covariates
at 1000 observations, seven of them at zero, it read 2.47 log-likelihood
units where the free coordinates' own reading is 4.5e-08, and a marginal
search whose every point was refused a resolution for that reason ran
out of backtracks and reported failure. Whether such a coordinate should
leave zero is the kinked block's own condition, read by
[`alternation_readings()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md).
