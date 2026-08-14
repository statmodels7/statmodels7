# One Setting of the Path, Read From the Term

One Setting of the Path, Read From the Term

## Usage

``` r
statmod_path_setting(spec, row, field, default, name = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- field:

  Which field of the penalty's entry to read.

- default:

  What the criterion asks for.

- name:

  The hyperparameter's name, or `NULL` where the setting is one per
  term.

## Value

A single number.
