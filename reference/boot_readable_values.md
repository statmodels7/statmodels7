# The Readable Quantities of Every Replicate

Returns the matrix whose columns are the quantities that
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
reports with `readable = TRUE`, evaluated at each replicate. It is the
matrix whose covariance is taken by
[`boot_readable_cov()`](https://statmodels7.github.io/statmodels7/reference/boot_readable_cov.md).

## Usage

``` r
boot_readable_values(object, spec, design, B)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  Its specification and design.

- B:

  The matrix of replicates, over the joint coordinates.

## Value

A matrix with a row for each replicate and the names of the readable
quantities, or `B` when every readable quantity is a coordinate.
