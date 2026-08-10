# Print One Block of a Summary

A header saying what the block is and what it spends, then its rows.

## Usage

``` r
print_block(b, digits = 4L)
```

## Arguments

- b:

  A block record, as
  [`summary_blocks`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
  returns.

- digits:

  Significant digits.

## Value

`NULL`, invisibly.

## Details

A row whose quantity was held fixed prints its value and blanks the
rest, rather than showing `NA` four times over: what the columns say is
that nothing estimated it, and the mark in the header says so once.
