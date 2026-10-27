# Degrees of Freedom Addressed by Coordinate Name

Degrees of Freedom Addressed by Coordinate Name

## Usage

``` r
sandwich_df_by_name(fn, names)

# S3 method for class 'statmod_df_fn'
print(x, ...)
```

## Arguments

- fn:

  The function of positions that
  [`cr2_df_closure()`](https://statmodels7.github.io/statmodels7/reference/cr2_df_closure.md)
  returns.

- names:

  The names of the coordinates that the sandwich covers, in the order of
  its columns.

- x:

  The function of coordinate names.

- ...:

  Unused.

## Value

A function of a character vector of names that returns the degrees of
freedom, `NA` for a name the variance does not cover. Its second
argument `contrast` is a matrix whose rows are named by coordinates and
whose columns are contrasts; the degrees of freedom of those contrasts
are returned, `NA` for a contrast that reads a coordinate the variance
does not cover.
