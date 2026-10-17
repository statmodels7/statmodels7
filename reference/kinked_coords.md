# Which Stacked Coordinates a Kinked Penalty Covers

Returns the positions of the coefficients whose penalty is not
differentiable everywhere – a lasso, a SCAD, an MCP – as an integer
vector.

## Usage

``` r
kinked_coords(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

An integer vector of stacked positions, possibly empty.

## Details

The function tests the coordinates in the way that
[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
tests the blocks to decide which of the two fitting schemes is used for
a block, and it reads the same predicate
([`penalty_has_kink()`](https://statmodels7.github.io/statmodels7/reference/penalty_has_kink.md)),
so that the two cannot disagree. A penalty that is twice differentiable
contributes nothing, so that a ridge, a random effect and a smooth are
testable.

The result follows from the enumeration of each penalty, so a partial
penalty (one covering some of a term's coordinates and not the rest)
names only what it covers, and a term is not classified as a whole. The
penalty of a structural term covers positions among the term's own
parameters and not among the design's, and contributes nothing here.

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
the enumeration that it reads, and
[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md),
which gives the same classification for a whole block.
