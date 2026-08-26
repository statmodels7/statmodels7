# How Many Values a Path Visits for One Hyperparameter

The grid size the term asked for, or the layer's fallback where it asked
for nothing.

## Usage

``` r
statmod_grid_size(spec, row, default)
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

A single integer.

## Details

How finely a hyperparameter is swept is a property of the term for the
same reason as whether it is swept at all: a penalized block of four
columns and one of four hundred want different grids, and a criterion
applies to every hyperparameter of the model at once, the smooth ones
included, and those are read at the mode instead of being swept, so a
criterion cannot know which kind it is looking at.
[`modelterms7::term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.html)
is where a term says so, and the value travels with the penalty's entry,
so one reached through a sub-term of a structural term carries it too.

## See also

[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
