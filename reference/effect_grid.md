# A Grid over the Variables of an Effect

Builds new data for an effect: `n_grid` points over the range of a
numeric variable, or its levels for a factor, crossed with the levels of
a factor among `cross` (when it has at most twelve), and every other
variable of `vars_all` held at its median or its first level.

## Usage

``` r
effect_grid(data, var, vars_all, n_grid, cross = vars_all)
```

## Arguments

- data:

  The data frame of the fit.

- var:

  The variable on the horizontal axis.

- vars_all:

  Every variable the term reads.

- n_grid:

  The number of points for a numeric variable.

- cross:

  The variables that may be crossed with the grid: those that enter the
  drawn contribution together with `var`.

## Value

A data frame with the attribute `"group"`, the name of the crossed
factor or `NULL`.
