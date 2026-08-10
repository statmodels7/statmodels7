# The Subspace a Marginal Criterion Integrates Over

`NULL` for REML, which integrates everything, and an orthonormal basis
of the penalty's range space for ML.

## Usage

``` r
integrated_basis(spec, design, kind)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- kind:

  `"reml"` or `"ml"`.

## Value

`NULL`, or a matrix with one column per integrated direction.

## Details

The total penalty is block diagonal over the terms, each owning its own
columns, so its null space is the direct sum of the unpenalized columns
and each penalty's own null space. It is assembled that way and never
read off the eigenvalues of the assembled matrix, which measure the
arithmetic rather than the family: with two smoothing parameters ten
orders of magnitude apart, an eigenvalue count of a sum falls while the
null space does not move.
