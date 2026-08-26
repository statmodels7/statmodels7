# Record Where a Path Has Just Been

Writes the size of each kinked penalty's kink at the given
hyperparameters onto its block, so that the next point of a path can
screen against it.

## Usage

``` r
blocks_at_kink(blocks, hyper)
```

## Arguments

- blocks:

  The blocks, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns them, with a `sparse` list of the kinked entries.

- hyper:

  The hyperparameters at the point just fitted.

## Value

`blocks`, with each entry of its `sparse` list carrying a `prev_kink`
element: the size of that penalty's kink at `hyper`. The `smooth` list
is untouched.

## Details

The previous point travels on the blocks themselves, sparing every layer
between the path and the descent an argument. Where a penalty was a
moment ago is a property of its block, and the path rebuilds the blocks
at each point in any case.

## See also

[`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md),
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
