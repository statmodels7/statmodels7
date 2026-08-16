# The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated

\\1/\lVert A^{-1}\rVert_1\\ from a sparse Cholesky factor, which is the
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

The sparse route needs a condition estimate OF ITS OWN, and it cannot
borrow the dense one: `chol_rcond_cpp` reads a dense triangular factor.
[`Matrix::rcond`](https://rdrr.io/pkg/Matrix/man/rcond-methods.html) is
not the answer either – measured, it costs 10.3 ms at p = 503 and 500 ms
at p = 2003, more than the dense factorization the sparse route exists
to replace. Higham's one-norm estimator applied to the factor's own
solves costs 0.58 ms at p = 53 and 0.80 ms at p = 1003, nearly flat,
because it is a handful of triangular solves and an R loop around them.

For a symmetric matrix \\\lVert A^{-1}\rVert_1 \ge \lVert A^{-1}\rVert_2
= 1/\lambda\_{\min}\\, so the quantity returned is at or below the
smallest eigenvalue, and the estimator's own error is a further
underestimate of the norm in the other direction. It is used exactly as
the dense estimate is: to separate a matrix comfortably positive
definite from one that is not, two situations that differ by some
fifteen orders of magnitude here (measured on a design with two
identical columns, 1.4e-14 relative to the matrix's scale, against 1.5e3
for a hyperparameter driven to 1e15). A factor of two either way cannot
move that verdict, which is the argument already recorded for the dense
estimator.

## See also

[`pd_factor`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md),
[`pd_logdet`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
