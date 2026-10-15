# The Log-Determinant of a Penalized Information, Robustly

\\\log\|M\|\\ for a matrix that a Laplace approximation requires to be
positive definite, by the cheap route where that is safe and a more
expensive one where it is not.

## Usage

``` r
pd_logdet(M, scale = NULL)
```

## Arguments

- M:

  A symmetric matrix.

- scale:

  A reference magnitude: the unpenalized information's own scale, so
  that a hyperparameter legitimately sent to 1e15 is told apart from a
  flat direction. `NULL` uses the 1-norm of `M`, and a larger value is
  capped at that norm.

## Value

A list with `logdet` and `ok`, or `ok = FALSE` where the matrix is not
positive definite.

## Details

**Why not [`chol()`](https://rdrr.io/r/base/chol.html) alone.** At a
condition number near the rounding floor, whether the factorization
raises an error is decided by the arithmetic and not by the matrix.
Reading the determinant off `chol(M)` alone would therefore report the
criterion as nonexistent at points where it exists, and the outer search
would backtrack from them. The same consideration applies to
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
and to the rank tests of basis7.

**The three routes.** The factorization is tried first, being O(p^3/3)
and the common case. Where it succeeds, LAPACK's condition estimator
reads the smallest eigenvalue off the factor already in hand for O(p^2),
and a matrix comfortably away from the floor is accepted with the
determinant that the factor gives. The eigendecomposition is computed
only where that test is inconclusive, or where the factorization raised
an error. It costs more and answers about the matrix: a factorization
that failed by rounding on a matrix that is in fact positive definite is
recovered there, and a matrix that is genuinely rank deficient is
rejected deterministically and not according to the platform.

## See also

[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
