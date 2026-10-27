# Whether a Test Refits the Model Under the Restriction

Returns `TRUE` for a statistic that is evaluated at the estimates in
which the tested coefficient is held at a value
([`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
and
[`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)),
and `FALSE` for the statistic that reads the unrestricted fit alone
([`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md)).

## Usage

``` r
test_restricted(test, ...)
```

## Arguments

- test:

  A
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
  object.

- ...:

  Unused.

## Value

A single logical.

## Details

A restricted statistic is about one coefficient, since one coordinate is
held.
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
reads this generic to require `readable = FALSE` for such a statistic,
and the summary reads it to decide whether a row costs a refit.

## See also

[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the four statistics and what each needs.

## Examples

``` r
test_restricted(wald_test())
#> [1] FALSE
test_restricted(lr_test())
#> [1] TRUE
```
