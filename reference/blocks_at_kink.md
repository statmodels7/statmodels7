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

`blocks`, in which each entry of its `sparse` list carries a `prev_kink`
element: the size of that penalty's kink at `hyper`. An entry whose kink
size is not finite and positive, or cannot be read, carries no
`prev_kink` element. The `smooth` list is untouched.

## Details

The size at the previous point travels on the blocks themselves, so that
the layers between the path and the descent do not have to pass an
argument.

## See also

[`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md),
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
