# The Working Weights of an Equation

Returns the information of the log-likelihood in the predictor of an
equation, \\w_i\\, expected or observed according to `information`, and
the working residuals \\s_i / w_i\\ of
[`working_residuals()`](https://statmodels7.github.io/statmodels7/reference/working_residuals.md),
from one evaluation of the derivatives.

## Usage

``` r
working_parts(fit, param, information)
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

A list with `wr` and `w`, numeric vectors with one value per
observation; `wr` is `NA` and `w` is zero where the information is not
positive.
