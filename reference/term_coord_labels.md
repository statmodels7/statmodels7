# The Coordinates of a Term's Own Covariance

The within-group column names of a random-effect term, which are what
the coordinates of its own multivariate prior are, or nothing for a term
that has no grouping.

## Usage

``` r
term_coord_labels(term, ii = NULL)
```

## Arguments

- term:

  One built term.

- ii:

  The entry's positions in the term's block, or `NULL` for every
  within-group column.

## Value

A character vector, possibly empty.

## Details

Where the prior covers some of the within-group columns and not all of
them, as the null-space prior of a smooth random effect does, `ii` gives
the entry's positions in the block and the labels are those of the
columns it covers, in the order of the block.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
for the same question about a shared block,
[`entry_owner()`](https://statmodels7.github.io/statmodels7/reference/entry_owner.md)
for the term to which this applies.
