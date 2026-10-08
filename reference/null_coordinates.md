# The Null Coordinates of a Penalty

The positions, within a penalty's block, that span its null space, or
`NULL` where the null space is not spanned by coordinates.

## Usage

``` r
null_coordinates(pen, k)
```

## Arguments

- pen:

  A penalties7 penalty.

- k:

  The block's width.

## Value

An integer vector of positions within the block, empty for a penalty of
full rank, or `NULL`.
