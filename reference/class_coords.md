# What Each Coordinate of a Covariance Block Is

One row per coordinate of a class's covariance, naming the equation it
belongs to, the term that carries it, the within-group column it is, and
the label a summary prints for it.

## Usage

``` r
class_coords(u)

coord_labels_of(cd)
```

## Arguments

- u:

  One penalized unit carrying `pieces`, as
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  returns for a covariance class.

- cd:

  The first three columns of the result.

## Value

A data frame with `param`, `term`, `column` and `label`, one row per
coordinate.

## Details

The prior of a covariance class describes the effects of **one group**
over every column the label collects, so its dimension is the sum of the
members' within-group widths and its coordinates run over the members in
order, each contributing its own columns in order – which is exactly how
[`class_index()`](https://statmodels7.github.io/statmodels7/reference/class_pieces.md)
interleaves them.

The multivariate family that carries the chart numbers those coordinates
`v1`, `v2`, and so on, and it is right to: it is a law on
\\\mathbb{R}^d\\ and knows nothing of the model. Which equation and
which term a coordinate belongs to is this layer's answer, and without
it a printed correlation between `v1` and `v3` says nothing at all.

The label is the shortest one that separates the coordinates: the column
alone where every coordinate comes from one term, the equation and the
column where more than one term is involved, and the term as well where
two terms of one equation write the same column name. A coordinate
reached through a subformula carries the parameter it develops in front
of its column, so that an effect on a break-point is not read as an
effect on the equation the break-point sits in.

## See also

[`class_notes()`](https://statmodels7.github.io/statmodels7/reference/class_notes.md),
[`summary_class_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_class_blocks.md),
which print it;
[`modelterms7::term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.html),
which supplies the column names.
