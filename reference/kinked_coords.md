# Which Stacked Coordinates a Kinked Penalty Covers

The positions of the coefficients whose penalty is not differentiable
everywhere – a lasso, a SCAD, an MCP – as an integer vector.

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

It is the question
[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
asks to decide which of the two fitting schemes a block goes to, asked
of the coordinates instead of the blocks, and it reads the same
predicate
([`penalty_has_kink()`](https://statmodels7.github.io/statmodels7/reference/penalty_has_kink.md))
so the two cannot disagree. A penalty that is twice differentiable
contributes nothing, which is what makes a ridge, a random effect and a
smooth testable.

A penalty's own enumeration is what answers, so a partial penalty – one
covering some of a term's coordinates and not the rest – names only what
it covers, and a term is not classified as a whole. A structural term's
penalty covers positions among the TERM'S own parameters rather than the
design's, and contributes nothing here.

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
the enumeration it reads, and
[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md),
which asks the same question of a whole block.
