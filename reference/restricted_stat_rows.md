# A Table's Statistics From a Restricted Fit

The likelihood-ratio, score or gradient statistic for each coefficient a
restricted fit can hold, against the null that it is zero, reported as
the signed root of the \\\chi^2_1\\ value together with that value's
p-value.

## Usage

``` r
restricted_stat_rows(fit, ci, test, type, spec, design)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ci:

  The flat interval table, with the Wald statistic already in it.

- test:

  Which statistic: `"lr"`, `"score"` or `"gradient"`.

- type:

  Which variance the Wald bracket reads, passed on.

- spec, design:

  The specification and its design.

## Value

A list with `statistic`, `p_value` and `note`, the last a character
vector of what the reader has to be told about the column.

## Details

The root carries the sign of the estimate, so the column holds one kind
of quantity whatever produced it and Wald's own \\z\\ is the case where
the root is exact. The p-value is the \\\chi^2_1\\ one and is not
recomputed from the root.

A row
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
does not name keeps `NA` in both columns rather than the Wald reading it
arrived with: a column carrying two different tests, one row apiece and
unmarked, is worse than a column with holes in it.

## See also

[`statmod_stat()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat.md),
which computes one of them.
