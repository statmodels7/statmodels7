# Which Coefficients a Search Should Cover

The positions, in the stacked coefficient vector, of the coordinates
where the likelihood is not convex: each equation's intercept, and the
blocks of the terms that recompute their own design.

## Usage

``` r
search_coords(spec, design)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

## Value

An integer vector of positions, possibly empty.

## Details

The predicate for the second is
[`refreshes_own_block`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md),
the same one
[`unfittable_reason`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md)
uses, so a term written later is covered without an edit here. A convex
block is left out deliberately rather than forgotten: the scoring step
reaches its optimum from anywhere, and a search over it would spend the
budget where it buys nothing.
