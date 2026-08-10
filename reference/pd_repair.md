# Floor the Eigenvalues of a Curvature Matrix

Returns the matrix unchanged when it is positive definite, and otherwise
its eigendecomposition with the eigenvalues floored.

## Usage

``` r
pd_repair(A, rel = 1e-08)
```

## Arguments

- A:

  A symmetric matrix.

- rel:

  The floor, relative to the largest eigenvalue.

## Value

A symmetric positive definite matrix.

## Details

Abandoning a start because the curvature is indefinite is what
[`solve()`](https://rdrr.io/r/base/solve.html) forces and what a repair
avoids; the floor is relative to the largest eigenvalue, since an
absolute one means nothing across scales.
