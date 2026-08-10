# A Square-Root Factor of the Penalty

Returns \\C\\ with \\C'C = S(\theta)\\, the block-diagonal penalty
Hessian, dropping the rows a null space contributes nothing to.

## Usage

``` r
penalty_sqrt(S)
```

## Arguments

- S:

  The penalty Hessian.

## Value

A matrix with `ncol(S)` columns, or `NULL`.

## Details

A penalty is positive semidefinite and may be rank deficient – a
spline's is, by exactly the dimension of its null space – so the factor
comes from an eigendecomposition with the non-positive eigenvalues
dropped, rather than from a Cholesky, which would fail there. A
non-convex penalty, whose Hessian is indefinite, has no such factor and
the caller falls back.
