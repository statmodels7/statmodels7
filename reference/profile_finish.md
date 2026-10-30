# The Legend of a Profile Panel

Adds the scales that name every element of a panel: the colour scale for
the drawn series, with the key of each entry drawn as its line or its
point, the line type of the estimate and of the limit of the interval,
and the shape of the points whose refit did not converge.

## Usage

``` r
profile_finish(p, entries)
```

## Arguments

- p:

  A ggplot2 object.

- entries:

  The entries of the colour scale, as
  [`profile_base()`](https://statmodels7.github.io/statmodels7/reference/profile_base.md)
  returns them.

## Value

`p` with the scales and the legend at the bottom.
