# S7 Class of the Influence Graph

The class that
[`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotInfluence(
  label = character(0),
  param = character(0),
  measure = character(0),
  adjust = logical(0),
  information = NULL
)
```

## Arguments

- label, param, measure, adjust, information:

  See
  [`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md).

## Value

A `PlotInfluence` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
