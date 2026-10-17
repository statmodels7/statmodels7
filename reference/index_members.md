# The Hyperparameters Behind Each Row of the Outer Index

Returns the member table: one line per estimated hyperparameter, with
the row of the index that stands for it.

## Usage

``` r
index_members(idx)
```

## Arguments

- idx:

  The index, as
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
  returns it.

## Value

A data frame with columns `row`, `parameter`, `term` and `name`.

## Details

If no hyperparameter is shared, the member table is the index itself
with a `row` column counting up. Where a label ties several
hyperparameters, they carry the same `row`, and a quantity accumulated
over members lands on the coordinate of the group.
