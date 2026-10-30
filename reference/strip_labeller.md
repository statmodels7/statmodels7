# The Headings of the Panels of a Graph

Returns a labeller that heads each panel with the name of the variable
that divides the panels and the value of the panel, as in "x:
(0.29,0.5\]".

## Usage

``` r
strip_labeller(name)
```

## Arguments

- name:

  The name of the variable, or `NULL` for the value alone.

## Value

A labeller for
[`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
