# S7 Class of the Calibration Graph

The class that
[`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotCalibration(
  label = character(0),
  probs = integer(0),
  over = NULL,
  n_inter = integer(0),
  level = integer(0),
  seed = NULL
)
```

## Arguments

- label, probs, over, n_inter, level, seed:

  See
  [`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md).

## Value

A `PlotCalibration` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
