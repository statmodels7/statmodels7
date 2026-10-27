# The Data of a Quantile-Quantile or Worm Graph

Sorts the quantile residuals (within the classes of `over`, when given)
and pairs them with the normal plotting positions and the pointwise
band; for a worm plot the plotting position is subtracted from the
residuals and from the band.

## Usage

``` r
diag_qq_data(type, fit, ctx, worm, over = NULL, n_inter = 4)
```

## Arguments

- type:

  A `PlotQq` or `PlotWorm` object.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- worm:

  Whether to detrend.

- over:

  The name of the covariate that divides the observations, or `NULL`.

- n_inter:

  The number of classes.

## Value

A list with the data frames `points` and `band`.
