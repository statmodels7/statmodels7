# S7 Class of the Break-Point Graph

The class that
[`plot_breakpoint()`](https://statmodels7.github.io/statmodels7/reference/plot_breakpoint.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotBreakpoint(
  label = character(0),
  profile = logical(0),
  test = NULL,
  level = integer(0),
  vcov_type = NULL,
  information = NULL,
  n_grid = integer(0),
  n_profile = integer(0)
)
```

## Arguments

- label, profile, test, level, vcov_type, information, n_grid,
  n_profile:

  See
  [`plot_breakpoint()`](https://statmodels7.github.io/statmodels7/reference/plot_breakpoint.md).

## Value

A `PlotBreakpoint` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
