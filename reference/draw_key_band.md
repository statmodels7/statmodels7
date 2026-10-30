# The Legend Key of an Estimate with Its Band

Draws the key of a line together with a light band behind it, in the
colour of the line, so that one legend entry names an estimate and its
interval. It is passed as `key_glyph` to the layer of the line, and the
layer of the band draws no key.

## Usage

``` r
draw_key_band(data, params, size)
```

## Arguments

- data, params, size:

  As for the `draw_key_*()` functions of ggplot2.

## Value

A grob.
