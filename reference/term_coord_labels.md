# The Coordinates of a Term's Own Covariance

The within-group column names of a random-effect term, which are what
the coordinates of its own multivariate prior are, or nothing for a term
that has no grouping.

## Usage

``` r
term_coord_labels(term)
```

## Arguments

- term:

  One built term.

## Value

A character vector, possibly empty.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
for the same question about a shared block,
[`entry_owner()`](https://statmodels7.github.io/statmodels7/reference/entry_owner.md)
for the term this is asked of.
