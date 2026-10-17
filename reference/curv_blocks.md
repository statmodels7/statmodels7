# The Blocks With the Penalties That a Specification Now Carries

Replaces the penalty of every kinked block with the one that the
specification carries under the same key, leaving everything else that
the block records where it was.

## Usage

``` r
curv_blocks(blocks, spec, design)
```

## Arguments

- blocks:

  The blocks, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns them.

- spec, design:

  The specification carrying the new penalties, and its design.

## Value

`blocks`, in which the `penalty` of each kinked entry is replaced.

## Details

[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
rewrites a scaled SCAD or MCP's curvature at every pass. Rebuilding the
blocks with
[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
would lose what a path has written onto them, including the previous
point's kink, so only the penalty is replaced.

## See also

[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
