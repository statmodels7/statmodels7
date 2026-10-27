# S7 Class of the Quantile-Quantile Graph

The class that
[`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotQq(
  label = character(0),
  band = character(0),
  level = integer(0),
  n_sim = integer(0),
  seed = NULL
)
```

## Arguments

- label, band, level, n_sim, seed:

  See
  [`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md).

## Value

A `PlotQq` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
