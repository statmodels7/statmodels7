# The Values a Term Wrote Out for One Hyperparameter

The grid the TERM named, or `NULL` where it left the path to build one.

## Usage

``` r
statmod_values(spec, row)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

## Value

A numeric vector, or `NULL`.

## Details

The third state of a hyperparameter's argument, beside holding it at one
number and leaving it to be estimated over a grid the path constructs.
It travels with the penalty's entry like the grid size and the depth, so
a penalty reached through a sub-term of a structural term carries it
too.

## See also

[`term_values`](https://statmodels7.github.io/modelterms7/reference/term_values.html),
[`path_forced`](https://statmodels7.github.io/statmodels7/reference/path_forced.md)
