# S7 Class of the Surface Graph

The class that
[`plot_surface()`](https://statmodels7.github.io/statmodels7/reference/plot_surface.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotSurface(
  label = character(0),
  legend = TRUE,
  mask = logical(0),
  slices = NULL,
  level = integer(0),
  vcov_type = NULL,
  n_grid = integer(0)
)
```

## Arguments

- label, mask, slices, level, vcov_type, n_grid:

  See
  [`plot_surface()`](https://statmodels7.github.io/statmodels7/reference/plot_surface.md).

## Value

A `PlotSurface` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
