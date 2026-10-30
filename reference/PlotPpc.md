# S7 Class of the Posterior Predictive Check

The class that
[`plot_ppc()`](https://statmodels7.github.io/statmodels7/reference/plot_ppc.md)
instantiates. Use the constructor, which checks the properties.

## Usage

``` r
PlotPpc(
  label = character(0),
  legend = TRUE,
  statistic = NULL,
  n_sim = integer(0),
  seed = NULL
)
```

## Arguments

- label, statistic, n_sim, seed:

  See
  [`plot_ppc()`](https://statmodels7.github.io/statmodels7/reference/plot_ppc.md).

## Value

A `PlotPpc` object, inheriting from
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
