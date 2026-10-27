# The Profile Inset of the Break-Point Graph

Builds the layers of the inset: a small ggplot of the profile, placed
with
[`ggplot2::annotation_custom()`](https://ggplot2.tidyverse.org/reference/annotation_custom.html)
in a band above the data, and a blank layer that extends the vertical
axis to make room for it. The inset is on the side of the graph, left or
right, where the data above their middle are fewer.

## Usage

``` r
breakpoint_inset(data)
```

## Arguments

- data:

  The list returned by the
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
  method of
  [`plot_breakpoint()`](https://statmodels7.github.io/statmodels7/reference/plot_breakpoint.md).

## Value

A list of ggplot2 components.
