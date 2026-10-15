# Floor the Eigenvalues of a Curvature Matrix

Returns the matrix unchanged when it is positive definite, and otherwise
the matrix rebuilt from its eigendecomposition with the eigenvalues
floored.

## Usage

``` r
pd_repair(A, rel = 1e-08)
```

## Arguments

- A:

  A symmetric matrix.

- rel:

  The floor, relative to the largest absolute eigenvalue.

## Value

A symmetric positive definite matrix, or `NULL` where `A` has a
non-finite entry. A matrix with no rows is returned as it is.

## Details

Abandoning a start because the curvature is indefinite is what
[`solve()`](https://rdrr.io/r/base/solve.html) forces and what a repair
avoids; the floor is relative to the largest absolute eigenvalue, since
an absolute one means nothing across scales. Positive definiteness is
tested by a Cholesky factorization.
