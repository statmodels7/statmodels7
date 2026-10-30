# S7 Class of the Derivative Graph

The class that
[`plot_derivative()`](https://statmodels7.github.io/statmodels7/reference/plot_derivative.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotDerivative(
  label = character(0),
  legend = TRUE,
  order = integer(0),
  level = integer(0),
  simultaneous = logical(0),
  n_sim = integer(0),
  seed = NULL,
  vcov_type = NULL,
  n_grid = integer(0)
)
```

## Arguments

- label, order, level, simultaneous, n_sim, seed, vcov_type, n_grid:

  See
  [`plot_derivative()`](https://statmodels7.github.io/statmodels7/reference/plot_derivative.md).

## Value

A `PlotDerivative` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
