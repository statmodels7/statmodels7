# S7 Class of the Effect Graph

The class that
[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotEffect(
  label = character(0),
  residuals = logical(0),
  binned = NULL,
  scale = character(0),
  at = NULL,
  rug = logical(0),
  level = integer(0),
  vcov_type = NULL,
  information = NULL,
  n_grid = integer(0)
)
```

## Arguments

- label, residuals, binned, scale, at, rug, level, vcov_type,
  information, n_grid:

  See
  [`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md).

## Value

A `PlotEffect` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
