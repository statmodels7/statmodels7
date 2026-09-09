# The Stem Shared by a Block's Coefficient Names

Returns the part every coefficient name of a block has in common before
its own suffix, which is what identifies the term across two models
written on different bases.

## Usage

``` r
block_stem(nms)
```

## Arguments

- nms:

  A block's coefficient names.

## Value

A single string, or `NA_character_` where the block has no stem.

## Details

A block is keyed in the design by the term's deparsed call, so
`s(x, bspline_smooth(k = 6))` and `s(x, bspline_smooth(k = 10))` are two
different keys and never meet by name. Their coefficients, on the other
hand, are `s(x).lin`, `s(x).z1`, ... on both sides, so the stem `s(x)`
names the term without naming the basis. The suffix is taken from the
LAST dot rather than the first, so a covariate whose own name carries
one – `s(my.var)` – keeps it.

`NA` is returned where the names carry no suffix at all, which is the
parametric block, whose columns are variables and levels and which is
matched column by column instead.

## See also

[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
whose method pairs blocks with it.
