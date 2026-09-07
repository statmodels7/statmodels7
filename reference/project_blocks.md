# Project a Reference Fit's Predictor onto the Blocks That Did Not Match

Estimates the coefficients of the blocks
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)'s
`StartFrom` method left pending – a term the two models share written on
a different basis – by least squares against the part of the reference's
predictor the blocks already carried across do not explain.

## Usage

``` r
project_blocks(spec, design, ref, rdesign, out, pending)
```

## Arguments

- spec:

  The
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  being started.

- design:

  Its design.

- ref:

  The reference
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- rdesign:

  The reference's design.

- out:

  The starting coefficients as they stand, with the matched blocks
  already written in.

- pending:

  A list, one entry per distribution parameter, of `list(nm, idx)` pairs
  naming a block and the columns it occupies.

## Value

A list with `start`, the coefficients with the pending blocks estimated,
and `taken`, a list of one-row-per-coefficient data frames recording
what was projected.

## Details

Two bases of the same term span different subspaces, so
\\B\_{\mathrm{new}}\beta\_{\mathrm{new}} = B\_{\mathrm{old}}
\beta\_{\mathrm{old}}\\ has no solution in general and the projection
\$\$\hat\beta_S = \arg\min\_\beta \lVert X_S\beta - r\rVert^2, \qquad r
= X^{\mathrm{ref}}\beta^{\mathrm{ref}} - X\_{-S}\beta\_{-S},\$\$ is what
remains: the best the new basis can do at reproducing the fitted
function. Only the pending columns \\S\\ are estimated, and \\r\\ is the
partial residual of the reference's predictor after the blocks already
matched by name are removed at the values they were given. Where a block
matched exactly there is nothing to solve, so a refit whose bases all
agree performs no arithmetic here and returns the coefficients it was
handed.

The projection is exact when the reference's function lies in the span
of the new basis, and reduces to the identity when the two bases
coincide, which is why the name match above is a fast path rather than a
different answer. Whatever the reference explains and the new model has
no column for – a covariate the formula dropped – stays in \\r\\ and is
absorbed by the pending blocks, which is the closest the new design can
come to the predictor it is being started from.

Three conditions gate it, each returning the coefficients unchanged. The
two responses must be identical, which is what says the two designs are
read at the same rows: the projection compares two predictors pointwise
and means nothing across different data. Neither model may carry a
structural term, whose contribution is a recursion's state rather than
\\X\beta\\. And neither design may carry an adjustment for that
parameter, which is the same condition for a block that moves with its
coefficients.

## See also

[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
whose method calls it.
