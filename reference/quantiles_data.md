# The Data of a Quantile Graph

Computes the quantile curves and the observed points of
[`plot_quantiles()`](https://statmodels7.github.io/statmodels7/reference/plot_quantiles.md).

## Usage

``` r
quantiles_data(type, x, fit, ctx, discrete)
```

## Arguments

- type, x, fit, ctx, discrete:

  The graph type, the family, the fit, the context of
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
  and whether the family is discrete.

## Value

A list with the data frames `curves` and `points`.
