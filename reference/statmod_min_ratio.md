# How Far Down the Path Reaches for One Term

The depth the TERM asked for, or the criterion's default where it asked
for nothing.

## Usage

``` r
statmod_min_ratio(spec, row, default)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- default:

  The criterion's own ratio.

## Value

A single number.

## Details

One number per term rather than one per hyperparameter, because only the
sweep by kink size uses it: a bounded hyperparameter is swept over its
own interval and a shape that does not move the kink over a geometric
grid above its lower bound.

## See also

[`statmod_grid_size`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md),
[`path_values`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
