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

A block is matched by the name the formula gave its term: column by
column for the parametric block, and as a whole for every other where
the coefficient names agree. What is left over – the same term written
on another basis, paired by
[`block_stem()`](https://statmodels7.github.io/statmodels7/reference/block_stem.md)
and
[`same_term_kind()`](https://statmodels7.github.io/statmodels7/reference/same_term_kind.md)
– goes to
[`project_blocks()`](https://statmodels7.github.io/statmodels7/reference/project_blocks.md),
which estimates it against the partial residual of the reference's
predictor. See
[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md)
for why the three rules differ.

The result carries an attribute `"taken"`, a data frame naming every
coefficient the reference answered for and whether it was `matched` or
`projected`, so a caller can see what was reused rather than infer it
from the fit that follows.
