# The Labelled Points of a Diagnostic Graph

Draws the points returned by
[`diag_top()`](https://statmodels7.github.io/statmodels7/reference/diag_top.md)
in the reference colour, with their row numbers above them, and the
legend that names the element drawn for every observation.

## Usage

``` r
diag_top_layers(top, all, all_colour)
```

## Arguments

- top:

  The data frame returned by
  [`diag_top()`](https://statmodels7.github.io/statmodels7/reference/diag_top.md).

- all:

  The name of the element drawn for every observation.

- all_colour:

  Its colour.

## Value

A list of ggplot2 components.
