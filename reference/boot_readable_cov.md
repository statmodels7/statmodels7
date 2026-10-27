# The Covariance of the Readable Quantities of the Replicates

Evaluates the readable quantities at every replicate and returns their
covariance.

## Usage

``` r
boot_readable_cov(object, spec, design, B)
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

A square matrix with the names of
[`readable_vcov()`](https://statmodels7.github.io/statmodels7/reference/readable_vcov.md).

## Details

When every readable quantity is one coefficient, the replicates are
reordered. Otherwise
[`readable_joint()`](https://statmodels7.github.io/statmodels7/reference/readable_joint.md)
is evaluated at each replicate.
