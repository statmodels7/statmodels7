# Starting Values From Another Fit

The reference fit's estimates where the two models share a coefficient,
the fallback strategy's answer everywhere else.

## Arguments

- strategy:

  A `StartFrom` object.

- spec, design, obj, ...:

  As in
  [`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md).

## Value

A named list of numeric vectors, with the attribute `"taken"`.

## Details

A block is matched by the name that the formula gave its term: column by
column for the parametric block, and as a whole for every other block
where the coefficient names agree. The remaining blocks, which are the
same term written on another basis and paired by
[`block_stem()`](https://statmodels7.github.io/statmodels7/reference/block_stem.md)
and
[`same_term_kind()`](https://statmodels7.github.io/statmodels7/reference/same_term_kind.md),
go to
[`project_blocks()`](https://statmodels7.github.io/statmodels7/reference/project_blocks.md),
which estimates them against the partial residual of the reference's
predictor. See
[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md)
for the three rules.

The result carries an attribute `"taken"`, a data frame that names every
coefficient taken from the reference and records whether it was
`matched` or `projected`, so that the reused coefficients can be read
directly.
