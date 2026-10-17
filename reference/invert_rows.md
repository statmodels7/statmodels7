# Replace a Table's Limits by an Inverted Test

Runs
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md)
on each row for which a restricted fit is defined and writes the two
limits back, leaving the rest of the table as it was.

## Usage

``` r
invert_rows(fit, out, level, test, type)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- out:

  The table that
  [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  has built, already subset.

- level:

  The confidence level.

- test:

  Which test to invert.

- type:

  Which variance sets the starting bracket.

## Value

`out`, in which the `lower` and `upper` columns are replaced.

## Details

A row that the restricted fit cannot hold keeps `NA` limits rather than
the Wald limits that it carried on arrival, so that one column never
carries two different intervals.
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
determines which rows those are; if it names none of the rows of the
table, the function signals an error, since a table that consists only
of missing values is less informative than the reason.

## See also

[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which does the search.
