# One-Step Deletion Quantities of an Equation

Computes, for each observation, the first-order change of the
coefficients of an equation on deleting it, its leverage and Cook's
distance, as described at
[`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md).

## Usage

``` r
influence_deltas(fit, param, information, adjust = TRUE)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The equation.

- information:

  An object from
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  or
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html).

- adjust:

  Whether to divide the step by one minus the leverage.

## Value

A list with `delta` (a matrix with one row per observation and one
column per coefficient of the equation that the variance reports),
`leverage`, `cook` and `p_eff`.
