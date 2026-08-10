# Where Each Stacked Coefficient Comes From

One row per stacked coefficient, naming the distribution parameter, the
term and the coefficient, and saying whether the term carries a penalty
and whether that penalty has a kink.

## Usage

``` r
coef_labels(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A data frame with as many rows as there are coefficients.
