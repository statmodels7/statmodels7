# Replace a Table's Limits by an Inverted Test

Runs
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md)
on each row a restricted fit is defined for and writes the two limits
back, leaving the rest of the table as it was.

## Usage

``` r
invert_rows(fit, out, level, method, type)
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

- method:

  Which test to invert.

- type:

  Which variance sets the starting bracket.

## Value

`out` with its `lower` and `upper` columns replaced.

## Details

A row the restricted fit cannot hold keeps `NA` limits rather than the
Wald ones it arrived with, so that one column never carries two
different intervals. Which rows those are is
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)'s
answer, asked once; where it is none of them the request is refused, a
table of nothing but missing values being a worse answer than the reason
for it.

## See also

[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which does the search.
