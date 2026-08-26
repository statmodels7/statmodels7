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
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- kind:

  `"reml"` or `"ml"`.

## Value

`NULL` for `"reml"`, which integrates every direction and needs no
basis. For `"ml"`, a `p x r` matrix with orthonormal columns spanning
the penalty's range space, `p` the total coefficient count and `r` its
rank.

## Details

The total penalty is block diagonal over the terms, each owning its own
columns, so its null space is the direct sum of the unpenalized columns
and each penalty's own null space. The basis is assembled that way.

It is never read off the eigenvalues of the assembled matrix. Such a
count measures the arithmetic and not the family: with two smoothing
parameters ten orders of magnitude apart the small contributions sink
below any tolerance and the count falls, while the null space itself
does not move.

## See also

[`penalty_range_basis()`](https://statmodels7.github.io/statmodels7/reference/penalty_range_basis.md)
for one penalty's contribution,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
for the distinction this implements.
