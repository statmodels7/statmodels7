# S7 Class of the Autocorrelation Graphs

The class that
[`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md)
and
[`plot_pacf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md)
instantiate. Use the constructors, which check the properties.

## Usage

``` r
PlotAcf(
  label = character(0),
  legend = TRUE,
  lag_max = integer(0),
  squared = logical(0),
  by = NULL,
  time = NULL,
  seed = NULL,
  partial = logical(0)
)
```

## Arguments

- label, legend, lag_max, squared, by, time, seed:

  See
  [`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md).

- partial:

  Whether the graph shows the partial autocorrelation.

## Value

A `PlotAcf` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
