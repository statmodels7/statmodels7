# The Penalized Block, Kept on the Design

Returns the
[`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
of an equation's columns and keeps it in the design's `eta_memo`
environment, so that a descent called many times over one design copies
the block once.

## Usage

``` r
coord_block_at(design, p, X, cols)
```

## Arguments

- design:

  The design.

- p:

  The equation, a parameter name.

- X:

  The equation's design, `design[[p]]$X`.

- cols:

  The block's column positions within it.

## Value

What
[`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
returns for `X` and `cols`.

## Details

Only a design whose blocks do not move with the coefficients carries
that environment; for any other design, a fresh
[`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
is built at each call, as before. The entry is keyed on the equation and
the columns.

## See also

[`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md),
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
