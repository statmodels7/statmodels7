# The Series of the Observations of a Fit

Returns the positions of the observations of each series in time order:
from `by` and `time` when the graph gives them, from the ordering that a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
or
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term recorded otherwise, and one series in row order when there is
neither.

## Usage

``` r
acf_series(type, fit, ctx)
```

## Arguments

- type:

  A `PlotAcf` object.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

A list of integer vectors.
