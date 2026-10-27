# The Replicates of a Model With a Structural Term

Returns the matrix of bootstrap replicates over the joint coordinates of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
as
[`boot_vcov_matrix()`](https://statmodels7.github.io/statmodels7/reference/boot_vcov_matrix.md)
reads it, with the aliased coordinates missing.

## Usage

``` r
boot_replicate_matrix(object, vt)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  A `BootstrapVcov` object.

## Value

A matrix with one row for each replicate used.
