# A Table's Statistics From a Restricted Fit

Computes the likelihood-ratio, score or gradient statistic for each
coefficient that a restricted fit can hold, against the null that it is
zero, and reports it as the signed root of the \\\chi^2_1\\ value
together with the p-value of that value.

## Usage

``` r
restricted_stat_rows(fit, ci, test, vcov_type, spec, design)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ci:

  The flat interval table, with the Wald statistic already in it.

- test:

  A restricted
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md):
  [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
  [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
  or
  [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md).

- vcov_type:

  The
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  that the Wald bracket reads, passed on.

- spec, design:

  The specification and its design.

## Value

A list with `statistic`, `p_value` and `note`, the last a character
vector of the notes that qualify the column.

## Details

The root carries the sign of the estimate, so the column holds one kind
of quantity whichever statistic produced it, and Wald's own \\z\\ is the
case where the root is exact. The p-value is the \\\chi^2_1\\ one and is
not recomputed from the root.

A row that
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
does not name keeps `NA` in both columns instead of the Wald value with
which it arrived, so that each column carries a single test.

## See also

[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
which computes one of them.
