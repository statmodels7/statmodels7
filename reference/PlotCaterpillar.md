# S7 Classes of the Graphs of a Random Effect

The classes that
[`plot_caterpillar()`](https://statmodels7.github.io/statmodels7/reference/plot_caterpillar.md),
[`plot_prior()`](https://statmodels7.github.io/statmodels7/reference/plot_prior.md),
[`plot_random_cov()`](https://statmodels7.github.io/statmodels7/reference/plot_random_cov.md)
and
[`plot_curves()`](https://statmodels7.github.io/statmodels7/reference/plot_curves.md)
instantiate. Use the constructors, which check the properties.

## Usage

``` r
PlotCaterpillar(
  label = character(0),
  legend = TRUE,
  raw = logical(0),
  sort = logical(0),
  level = integer(0),
  vcov_type = NULL,
  information = NULL
)

PlotPrior(label = character(0), legend = TRUE)

PlotRandomCov(
  label = character(0),
  legend = TRUE,
  level = integer(0),
  which = integer(0)
)

PlotCurves(
  label = character(0),
  legend = TRUE,
  n_grid = integer(0),
  population = logical(0)
)
```

## Arguments

- label, raw, sort, level, vcov_type, information, which, n_grid,
  population:

  See the constructors.

## Value

An object inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
