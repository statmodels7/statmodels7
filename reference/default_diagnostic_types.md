# The Four Diagnostic Graphs Drawn by Default

Returns the graph types that
[`diagnostics()`](https://statmodels7.github.io/statmodels7/reference/diagnostics.md)
composes when `plot_types` is not given:
[`plot_ppc()`](https://statmodels7.github.io/statmodels7/reference/plot_ppc.md),
[`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md),
[`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md)
and
[`plot_residuals()`](https://statmodels7.github.io/statmodels7/reference/plot_residuals.md)
against the predicted median.

## Usage

``` r
default_diagnostic_types()
```

## Value

A named list of four
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md)
objects.
