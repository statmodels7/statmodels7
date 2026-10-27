# Build a Bootstrap Interval Object

Checks the arguments of
[`percentile_interval()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md),
[`basic_interval()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md)
and
[`bca_interval()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md)
and builds the object.

## Usage

``` r
boot_interval_new(class, name, n, resample, refit)
```

## Arguments

- class:

  The S7 class.

- name:

  The short name of the interval.

- n, resample, refit:

  As in
  [`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md).

## Value

An object that inherits from
[bootstrap_interval](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md).
