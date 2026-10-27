# The Bootstrap Variance That Supplies a Bootstrap Interval

Returns the resolved
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
object whose entry the interval reads, computing the replicates that are
missing.

## Usage

``` r
boot_interval_vcov(object, interval)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- interval:

  A
  [bootstrap_interval](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md)
  object.

## Value

A `BootstrapVcov` object with no `NULL` property.

## Details

A `NULL` number of replicates takes that of the single entry that
matches the other settings, and 1000 when none matches.
