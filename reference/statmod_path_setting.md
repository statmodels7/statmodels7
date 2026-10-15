# One Setting of the Path, Read From the Term

Reads one setting of the path from the penalty's entry in the term, and
returns `default` where the term does not set it.

## Usage

``` r
statmod_path_setting(spec, row, field, default, name = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- field:

  Which field of the penalty's entry to read.

- default:

  The fallback of
  [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
  for this setting, used where the term specifies none.

- name:

  The hyperparameter's name, or `NULL` where the setting is one per
  term.

## Value

A single number.
