# The Legend of the Elements of a Graph

Returns the scales that name the elements of a graph in a legend below
the panel. Each layer maps a constant string, the description of the
element, to `colour` or `fill`, and `colour` and `fill` are named
vectors whose names are those strings and whose values are the colours
drawn. A string that no layer maps does not appear in the legend; the
others appear in the order given.

## Usage

``` r
plot_key(colour = NULL, fill = NULL)
```

## Arguments

- colour, fill:

  Named character vectors, or `NULL`.

## Value

A list of ggplot2 components.
