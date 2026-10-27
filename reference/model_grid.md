# New Data on a Grid of One Covariate

Builds the grid of `over` (`n_grid` points over the range of a numeric
covariate, the levels of a factor) and completes it with the other
variables of the model, held at `at` or at their median or first level.

## Usage

``` r
model_grid(type, fit, ctx)
```

## Arguments

- type:

  A graph type with the property `over`.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

A list with `nd` (the new data) and `x` (the values of `over`).
