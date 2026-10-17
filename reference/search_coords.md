# Which Coefficients a Search Should Cover

Returns the positions, in the stacked coefficient vector, of the
coordinates where the likelihood is not convex: each equation's
intercept, and the blocks of the terms that recompute their own design.

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
[`refreshes_own_block()`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md),
which tests whether the term registers a
[`modelterms7::term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html)
method of its own, so a term written later is covered with no edit here.
A convex block is left out deliberately: the scoring step reaches its
optimum from any start, and a search over it would use the budget where
it gains nothing.

A penalized coordinate is then removed wherever it sits, which differs
from testing whether its term is convex: a penalty declared through a
sub-formula of a break-point or nonlinear term covers columns of the
block of that term, and the loop takes the whole block.
