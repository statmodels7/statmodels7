# The Coefficients an Axis Covers, Across Its Members

The union of the blocks the axis's members penalize, as positions in the
stacked coefficient vector.

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
block is empty has to be decided over all of them at once. Where nothing
is shared, the result is the single block's own index.
