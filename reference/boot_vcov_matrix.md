# The Covariance of the Stored Replicates

Returns the variance matrix that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
reports for
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md):
the covariance of the first `n` replicates of the entry, on the
coordinates of the other conventions.

## Usage

``` r
boot_vcov_matrix(object, vt, readable, param)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  A `BootstrapVcov` object.

- readable, param:

  As in
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A square matrix.

## Details

An aliased coordinate is missing in its row and column. A coefficient
under a kinked penalty that the fit holds at zero is included, with the
variance of its replicates. With `readable = TRUE` the quantities that
the model reports are computed for each replicate and their covariance
is taken: a term that reports a quantity other than its coefficients
([`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html))
is evaluated at every replicate, which is exact and not a delta method.
A readable quantity that reads a missing coordinate is missing.
