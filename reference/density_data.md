# The Data of a Density Graph

Computes the density or mass of each row of
[`plot_density()`](https://statmodels7.github.io/statmodels7/reference/plot_density.md).

## Usage

``` r
density_data(type, x, fit, ctx, discrete)
```

## Arguments

- type, x, fit, ctx, discrete:

  The graph type, the family, the fit, the context of
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
  and whether the family is discrete.

## Value

A list with the data frames `dens` and `points`.
