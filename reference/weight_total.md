# The Number of Observations Counted by the Weights

Returns the sum of the prior weights, which are frequency weights: a
weight of two counts the row twice. Without weights it is the number of
rows.

## Usage

``` r
weight_total(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A single number.

## See also

[`nobs.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/nobs.StatmodFit.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
