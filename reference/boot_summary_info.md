# What the Summary Reports About a Bootstrap Variance

Returns the number of replicates, how many were replaced, the relative
Monte Carlo error of a standard error in per cent, and the frequency
with which each coefficient under a kinked penalty was not zero.

## Usage

``` r
boot_summary_info(object, vt)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  A resolved `BootstrapVcov` object.

## Value

A list.
