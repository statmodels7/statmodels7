# S7 Classes of the Shrinkage and Path Graphs

The classes that
[`plot_shrinkage()`](https://statmodels7.github.io/statmodels7/reference/plot_shrinkage.md)
and
[`plot_path()`](https://statmodels7.github.io/statmodels7/reference/plot_path.md)
instantiate. Use the constructors, which check the properties.

## Usage

``` r
PlotShrinkage(label = character(0), sort = logical(0), information = NULL)

PlotPath(label = character(0), log_lambda = NULL, information = NULL)
```

## Arguments

- label, sort, information, log_lambda:

  See
  [`plot_shrinkage()`](https://statmodels7.github.io/statmodels7/reference/plot_shrinkage.md)
  and
  [`plot_path()`](https://statmodels7.github.io/statmodels7/reference/plot_path.md).

## Value

A `PlotShrinkage` or `PlotPath` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).
