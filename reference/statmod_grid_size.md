# How Many Values a Path Visits for One Hyperparameter

Returns the grid size that the term specifies, or the fallback of
[`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
where the term specifies none.

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

  The fallback of
  [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md),
  used for a term that specifies none.

## Value

A single integer.

## Details

How finely a hyperparameter is swept is a property of the term for the
same reason as whether it is swept at all: a penalized block of four
columns and one of four hundred require different grids. A criterion
applies to every hyperparameter of the model at once, including the
smooth ones, and those are read at the mode instead of being swept, so
the grid size cannot be a property of the criterion. A term gives its
grid size through
[`modelterms7::term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.html),
and the value travels with the penalty's entry, so one reached through a
sub-term of a structural term carries it too.

## See also

[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
