# Replace a Table's Limits by an Inverted Test

Runs
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md)
on each row a restricted fit is defined for and writes the two limits
back, leaving the rest of the table as it was.

## Usage

``` r
invert_rows(fit, out, level, test, type)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- out:

  The table
  [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  has built, already subset.

- level:

  The confidence level.

- test:

  Which test to invert.

- type:

  Which variance sets the starting bracket.

## Value

`out` with its `lower` and `upper` columns replaced.

## Details

A row the restricted fit cannot hold keeps `NA` limits rather than the
Wald ones it arrived with, so that one column never carries two
different intervals.
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
determines which rows those are; if it names none of the rows of the
table, the function signals an error, since a table of missing values
alone is less informative than the reason.

## See also

[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which does the search.
