# Quantile Residuals of a Response

Returns the quantile residuals of
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
for the response of the fit, or for another response `y` evaluated under
the fitted parameters (no refit).

## Usage

``` r
diag_resid(fit, y = NULL, seed = NULL)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- y:

  A numeric vector of the length of the response, or `NULL`.

- seed:

  Passed to
  [`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md).

## Value

A numeric vector.
