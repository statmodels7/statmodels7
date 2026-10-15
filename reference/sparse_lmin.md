# The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated

\\1/\lVert A^{-1}\rVert_1\\ from a sparse Cholesky factor, the same
quantity LAPACK's `dpocon` produces from a dense one.

## Usage

``` r
sparse_lmin(L, p)
```

## Arguments

- L:

  A `CHMfactor`.

- p:

  The order of the matrix.

## Value

A single number, or `NA_real_` where the estimate failed.

## Details

The sparse route needs a condition estimate of its own, because
`chol_rcond_cpp` reads a dense triangular factor.
[`Matrix::rcond`](https://rdrr.io/pkg/Matrix/man/rcond-methods.html) is
not used, since it costs more than the dense factorization that the
sparse route is meant to replace. The estimate is Higham's one-norm
estimator
([`Matrix::onenormest()`](https://rdrr.io/pkg/Matrix/man/condest.html))
applied to the solves with the factor, which costs a handful of
triangular solves.

For a symmetric matrix \\\lVert A^{-1}\rVert_1 \ge \lVert A^{-1}\rVert_2
= 1/\lambda\_{\min}\\, so the quantity returned is at or below the
smallest eigenvalue, and the error of the estimator is a further
underestimate of the norm in the other direction. It is used as the
dense estimate is: to separate a matrix comfortably positive definite
from one that is not, two situations whose smallest eigenvalues differ
by many orders of magnitude. An error of a factor of two in the estimate
does not change that verdict.

## See also

[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md),
[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
