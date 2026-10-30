# S7 Class of the Filter Graph

The class that
[`plot_filter()`](https://statmodels7.github.io/statmodels7/reference/plot_filter.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotFilter(
  label = character(0),
  legend = TRUE,
  level = integer(0),
  scale = character(0),
  response = logical(0),
  groups = NULL,
  vcov_type = NULL
)
```

## Arguments

- label, level, scale, response, groups, vcov_type:

  See
  [`plot_filter()`](https://statmodels7.github.io/statmodels7/reference/plot_filter.md).

## Value

A `PlotFilter` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
