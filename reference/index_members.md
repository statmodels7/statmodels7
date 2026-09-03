# The Hyperparameters Behind Each Row of the Outer Index

The member table: one line per estimated hyperparameter, with the row of
the index that stands for it.

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

Where nothing is shared it is the index itself with a `row` column
counting up, so a loop over members is the loop over rows that was there
before sharing existed. Where a label ties several, they carry the same
`row`, and a quantity accumulated over members lands on the group's
coordinate.
