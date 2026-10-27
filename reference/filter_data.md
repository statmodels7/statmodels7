# The Data of a Filter Graph

Computes the path of the filtered parameter of each equation in `params`
with its band, by group, and the observed response.

## Usage

``` r
filter_data(type, bp, params, fit, ctx)
```

## Arguments

- type:

  A `PlotFilter` object.

- bp:

  The blueprint of the term, or of its first component.

- params:

  The equations whose paths are drawn.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

A list with the data frames `path` and `obs`, and the attribute
`"labels"`.
