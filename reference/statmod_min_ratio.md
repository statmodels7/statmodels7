# How Far Down the Path Reaches for One Term

The depth the term asked for, or the layer's fallback where it asked for
nothing.

## Usage

``` r
statmod_min_ratio(spec, row, default)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- default:

  [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)'s,
  for a term that named none.

## Value

A single number.

## Details

One number per term, not one per hyperparameter, because only the sweep
by kink size reads it. A bounded hyperparameter is swept over its own
interval, and a shape that does not move the kink over a geometric grid
above its lower bound; neither takes its length from here.

## See also

[`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
