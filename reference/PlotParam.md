# S7 Classes of the Graphs of the Fitted Distribution

The classes that
[`plot_param()`](https://statmodels7.github.io/statmodels7/reference/plot_param.md),
[`plot_quantiles()`](https://statmodels7.github.io/statmodels7/reference/plot_quantiles.md),
[`plot_density()`](https://statmodels7.github.io/statmodels7/reference/plot_density.md)
and
[`plot_moments()`](https://statmodels7.github.io/statmodels7/reference/plot_moments.md)
instantiate. Use the constructors, which check the properties.

## Usage

``` r
PlotParam(
  label = character(0),
  param = character(0),
  over = character(0),
  at = NULL,
  level = integer(0),
  vcov_type = NULL,
  n_grid = integer(0)
)

PlotQuantiles(
  label = character(0),
  over = character(0),
  probs = integer(0),
  at = NULL,
  points = logical(0),
  n_grid = integer(0)
)

PlotDensity(
  label = character(0),
  at = NULL,
  over = NULL,
  n_grid = integer(0),
  scale = integer(0)
)

PlotMoments(
  label = character(0),
  over = character(0),
  moments = character(0),
  at = NULL,
  n_grid = integer(0)
)
```

## Arguments

- label, param, over, at, level, vcov_type, n_grid, probs, points,
  moments, scale:

  See the constructors.

## Value

An object inheriting from
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
