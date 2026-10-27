# Fitted Parameters on New Data

Calls
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for every parameter with the random effects set to zero, and recycles
each parameter to the number of rows.

## Usage

``` r
model_theta(fit, nd)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- nd:

  A data frame.

## Value

A named list of numeric vectors.
