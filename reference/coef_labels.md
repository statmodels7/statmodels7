# Where Each Stacked Coefficient Comes From

Returns one row per stacked coefficient. Each row names the distribution
parameter, the term and the coefficient, and records whether the term
carries a penalty and whether that penalty has a kink.

## Usage

``` r
coef_labels(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A data frame with as many rows as there are coefficients.
