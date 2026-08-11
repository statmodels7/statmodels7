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
  [`statmod_blocks`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns them.

- hyper:

  The hyperparameters at the point just fitted.

## Value

The blocks, each sparse entry carrying `prev_kink`.

## Details

The previous point travels on the blocks rather than through the
argument list of every layer between the path and the descent. It is a
property of the block – where its penalty was a moment ago – and the
path rebuilds the blocks at each point anyway.

## See also

[`coord_screen`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md),
[`statmod_path`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
