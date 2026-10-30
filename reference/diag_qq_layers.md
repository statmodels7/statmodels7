# The Layers of a Quantile-Quantile or Worm Graph

Draws the band, the reference line, the points and a legend that names
them and, when the data carry a `class` column, one facet for each
class, headed by the covariate and the class.

## Usage

``` r
diag_qq_layers(data, worm, type, over = NULL)
```

## Arguments

- data:

  The list returned by
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- worm:

  Whether the graph is detrended.

- type:

  The `PlotQq` or `PlotWorm` object, which gives the kind and the level
  of the band.

- over:

  The covariate that defines the classes, or `NULL`.

## Value

A list of ggplot2 components.
