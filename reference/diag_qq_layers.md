# The Layers of a Quantile-Quantile or Worm Graph

Draws the band, the reference line, the points and, when the data carry
a `class` column, one facet for each class.

## Usage

``` r
diag_qq_layers(data, worm)
```

## Arguments

- data:

  The list returned by
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- worm:

  Whether the graph is detrended.

## Value

A list of ggplot2 components.
