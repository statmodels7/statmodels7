# Project a Reference Fit's Predictor onto the Blocks That Did Not Match

Estimates the coefficients of the blocks that the `StartFrom` method of
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
left pending, which are terms that the two models share written on a
different basis. The estimate is by least squares against the part of
the reference's predictor that the blocks already carried across do not
explain.

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
\beta\_{\mathrm{old}}\\ has no solution in general, and the projection
\$\$\hat\beta_S = \arg\min\_\beta \lVert X_S\beta - r\rVert^2, \qquad r
= \eta^{\mathrm{ref}} - \eta\_{-S},\$\$ gives the closest reproduction
of the fitted function that the new basis allows. Only the pending
columns \\S\\ are estimated. The vector \\r\\ is the partial residual of
the reference's predictor \\\eta^{\mathrm{ref}}\\ after the blocks
already matched by name are removed at the values they were given
(\\\eta\_{-S}\\ is the new model's predictor with the pending columns at
zero). Both predictors are read through
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md),
so a block that moves with its coefficients (a break-point or a
nonlinear term) contributes its value and not its linearization
\\X\beta\\. Where every block matched exactly there is nothing to solve,
and a refit whose bases all agree returns the coefficients it was given.

The projection is exact when the reference's function lies in the span
of the new basis, and it reduces to the identity when the two bases
coincide, so the name match above is a shortcut and not a different
answer. Whatever the reference explains and the new model has no column
for, such as a covariate that the formula dropped, stays in \\r\\ and is
absorbed by the pending blocks.

Three conditions restrict the projection. The two responses must be
identical, so that the two designs are read at the same rows: the
projection compares two predictors observation by observation and has no
meaning across different data. Neither model may carry a structural
term, whose contribution is the state of a recursion and not a function
of the coefficients. Where either of these two conditions fails, the
coefficients are returned unchanged. A pending block that itself moves
with its coefficients (for example `seg(x)` refitted as
`seg(x, npsi = 2)`) is not projected, because its coefficients are not
those of a linear block, and it keeps the values given by the fallback
strategy. The pending columns of an equation are also left to the
fallback when their design has more than \\5 \times 10^7\\ entries.

## See also

[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
whose method calls it.
