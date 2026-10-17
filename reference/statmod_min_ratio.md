# The Depth of a Path for One Term

Returns the depth that the term specifies, or the fallback of
[`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
where the term specifies none.

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

  The fallback of
  [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md),
  used for a term that specifies none.

## Value

A single number.

## Details

The depth is one number per term and not one per hyperparameter, because
only the sweep by kink size reads it. A bounded hyperparameter is swept
over its own interval, and a shape that does not move the kink over a
grid whose distances above the smallest admissible shape are geometric;
in both sweeps the length does not come from this function.

## See also

[`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
