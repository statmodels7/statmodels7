# The Contribution of Several Parts of an Equation

Evaluates the sum of the contributions of `parts` on `newdata` (or on
the fitting rows), with the standard error that the joint variance of
their coefficients gives.

## Usage

``` r
parts_contribution(fit, param, parts, newdata, vcov_type)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The equation.

- parts:

  A list of parts, each with `name`, `term` and `cols`.

- newdata:

  A data frame, or `NULL`.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

## Value

A list with `fit` and `se`.
