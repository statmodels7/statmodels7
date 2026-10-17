# The Coefficients That an Axis Covers, Across Its Members

Returns the union of the blocks that the axis's members penalize, as
positions in the stacked coefficient vector.

## Usage

``` r
path_member_index(blocks, row)
```

## Arguments

- blocks:

  The block split.

- row:

  One row of
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

## Value

An integer vector of positions.

## Details

A shared axis is one value multiplying several penalties, so whether the
block is empty has to be decided over all of them at once. If the axis
has a single member, the result is that block's own index.
