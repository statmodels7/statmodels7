# Build an Autocorrelation Graph Type

Checks the arguments of
[`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md)
and
[`plot_pacf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md)
and returns the `PlotAcf` object.

## Usage

``` r
acf_type(label, partial, lag_max, squared, by, time, seed, legend)
```

## Arguments

- label:

  The name of the graph.

- partial:

  Whether the graph shows the partial autocorrelation.

- lag_max, squared, by, time, seed, legend:

  See
  [`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md).

## Value

A `PlotAcf` object.
