# S7 Class of the Worm Plot

The class that
[`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotWorm(
  label = character(0),
  over = NULL,
  n_inter = integer(0),
  band = character(0),
  level = integer(0),
  n_sim = integer(0),
  seed = NULL
)
```

## Arguments

- label, over, n_inter, band, level, n_sim, seed:

  See
  [`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md).

## Value

A `PlotWorm` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
