# The Data of a Criterion Graph

Computes the curve of the outer criterion against the logarithm of one
hyperparameter of a term, the fitted point and, for a refit, the
effective degrees of freedom of the term.

## Usage

``` r
criterion_plot_data(type, x, fit, ctx, ...)
```

## Arguments

- type:

  A `PlotCriterion` object.

- x:

  The term.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- ...:

  Unused.

## Value

A list with the data frames `curve` and `fitted`, and the attributes
`"axis"` and `"labels"`.
