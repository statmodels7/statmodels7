# S7 Class of the Autocorrelation Graph

The class that
[`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotAcf(
  label = character(0),
  lag_max = integer(0),
  squared = logical(0),
  seed = NULL
)
```

## Arguments

- label, lag_max, squared, seed:

  See
  [`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md).

## Value

A `PlotAcf` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
