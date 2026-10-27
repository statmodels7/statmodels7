# The Covariance of the Replicates of a Model With a Structural Term

Returns the covariance that
[`boot_vcov_matrix()`](https://statmodels7.github.io/statmodels7/reference/boot_vcov_matrix.md)
reports for a model whose replicates carry the free parameters of a
structural term after the coefficients.

## Usage

``` r
boot_vcov_structural(object, B, readable, spec, design, nm)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- B:

  The matrix of replicates, over the joint coordinates.

- readable:

  Whether to return the reported quantities.

- spec, design:

  The specification and design of the fit.

- nm:

  The names of the coefficient coordinates.

## Value

A square matrix.

## Details

The coordinates are the joint ones of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).
With `readable = TRUE` the reported quantities (a persistence, a
transition probability) are evaluated at every replicate, from the
replicate's coefficients and parameters, and their covariance is taken.
With `readable = FALSE` the matrix is ordered by equation as
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
orders it.
