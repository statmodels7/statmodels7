# The Data of an Effect Graph

Computes the curve (or the estimates of the levels), the band, the
partial residuals and their binned means for one or more parts of an
equation.

## Usage

``` r
effect_data(type, parts, fit, ctx)
```

## Arguments

- type:

  A `PlotEffect` object.

- parts:

  A list of parts, each with `name`, `term` and `cols`.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

A list with the data frames `curve`, `points`, `bins` and `rug`, and the
attributes `"labels"` and `"kind"` (`"numeric"` or `"factor"`).
