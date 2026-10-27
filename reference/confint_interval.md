# The Bootstrap Intervals of a Table of Coefficients

Computes the table of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
for an `interval` argument that is not `NULL`.

## Usage

``` r
confint_interval(
  object,
  parm,
  level,
  readable,
  test,
  interval,
  vcov_given,
  ...
)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- parm, level, readable:

  As in
  [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

- test:

  The `test` argument of
  [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

- interval:

  A
  [bootstrap_interval](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md)
  object.

- vcov_given:

  Whether the caller gave `vcov_type`.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

The data frame of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

## Details

The table is that of the Wald route with the bootstrap standard error,
which supplies the rows, the estimates and the standard errors; its
limits are replaced by the quantiles of the replicates. The variance
convention `vcov_type` and a test other than Wald's do not apply to a
bootstrap interval, and giving one is an error.
