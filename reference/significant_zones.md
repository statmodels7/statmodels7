# The Zones in Which a Band Excludes Zero

Returns the maximal runs of consecutive grid points at which `flag` is
`TRUE`, as intervals of the variable.

## Usage

``` r
significant_zones(x, flag)
```

## Arguments

- x:

  The grid, increasing.

- flag:

  A logical vector as long as `x`.

## Value

A data frame with `xmin` and `xmax`, one row for each run.
