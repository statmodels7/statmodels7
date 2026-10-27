# The Working Residuals of an Equation

Returns \\z_i - \eta_i = s_i / w_i\\, the difference between the working
response of an equation and its predictor, where \\s_i\\ is the score of
the log-likelihood with respect to the predictor and \\w_i\\ the
information, expected or observed, on the same scale. For a gaussian
mean with the identity link it is the ordinary residual.

## Usage

``` r
working_residuals(fit, param, information)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The equation.

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).

## Value

A numeric vector with one value per observation; `NA` where the
information is not positive.
