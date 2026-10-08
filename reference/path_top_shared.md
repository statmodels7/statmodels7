# The Top of a Path Over Every Member of an Axis

The size of the kink, in the axis's own block's units, at which every
member of the axis has all its coefficients at the kink.

## Usage

``` r
path_top_shared(obj, beta, blocks, row, hyper, b)
```

## Arguments

- obj:

  The stacked objective.

- beta:

  The coefficients the scores are read at.

- blocks:

  The block split.

- row:

  One row of
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- hyper:

  The hyperparameters.

- b:

  The axis's own block, `path_block(blocks, row)`.

## Value

A single number, a size of `b`'s kink.

## Details

A shared axis is one value multiplying several penalties, which may sit
in different equations – a lasso on the mean and a lasso on the scale
with one `id`. Each member's top is read from its own score, with every
member held at its kink, and carried onto the value of the
hyperparameter through its own penalty; the axis starts at the largest
of those values, written back as a size of the first member's kink,
which is the size the rest of the path speaks. Two blocks of one
equation sharing an `id` therefore start where a single block over both
starts, the score of a column not depending on how the columns were
grouped. Where nothing is shared this is
[`path_null_score()`](https://statmodels7.github.io/statmodels7/reference/path_null_score.md)
of the one block.

## See also

[`path_null_score()`](https://statmodels7.github.io/statmodels7/reference/path_null_score.md),
[`path_member_index()`](https://statmodels7.github.io/statmodels7/reference/path_member_index.md)
