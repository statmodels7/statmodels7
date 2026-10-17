# Floor the Eigenvalues of a Curvature Matrix

Returns the matrix unchanged when it is positive definite, and otherwise
returns the matrix rebuilt from its eigendecomposition after the
eigenvalues are floored.

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

A start at an indefinite curvature has to be abandoned if
[`solve()`](https://rdrr.io/r/base/solve.html) is used, and the repair
avoids this. The floor is relative to the largest absolute eigenvalue,
since an absolute floor has no meaning across scales. Positive
definiteness is tested by a Cholesky factorization.
