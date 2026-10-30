# S7 Class of the Residual Graph

The class that
[`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotResiduals(
  label = character(0),
  legend = TRUE,
  against = character(0),
  param = NULL,
  smooth = logical(0),
  seed = NULL
)
```

## Arguments

- label, against, param, smooth, seed:

  See
  [`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md).

## Value

A `PlotResiduals` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
