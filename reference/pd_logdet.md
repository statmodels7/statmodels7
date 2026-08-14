# The Log-Determinant of a Penalized Information, Robustly

\\\log\|M\|\\ for a matrix a Laplace approximation needs to be positive
definite, by the cheap route where that is safe and a costlier one where
it is not.

## Usage

``` r
pd_logdet(M, scale = NULL)
```

## Arguments

- M:

  A symmetric matrix.

- scale:

  A reference magnitude, as
  [`solve_pd`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  takes: the unpenalized information's own scale, so that a
  hyperparameter legitimately sent to 1e15 is told apart from a flat
  direction.

## Value

A list with `logdet` and `ok`, or `ok = FALSE` where the matrix is not
positive definite.

## Details

**Why not [`chol()`](https://rdrr.io/r/base/chol.html) alone.** A
marginal criterion read the determinant off `chol(M)` and reported the
criterion as NONEXISTENT whenever the factorization raised. At a
condition number near the rounding floor whether it raises is decided by
arithmetic and not by the matrix: measured on a hierarchical
score-driven panel, \\K+S\\ had a smallest eigenvalue of 4.3e-11 against
a condition number of 8.0e15, and the outer search then backtracked
through a dozen points reported unavailable towards one that had been
available a moment earlier. The same doubt this package already records
for
[`solve_pd`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
and for basis7's rank tests.

**The three routes.** The factorization is tried first, being O(p^3/3)
and the common case. Where it succeeds, LAPACK's condition estimator
reads the smallest eigenvalue off the factor already in hand for O(p^2),
and a matrix comfortably away from the floor is accepted with the
determinant the factor gives. Only where that test is inconclusive – or
the factorization raised at all – is the eigendecomposition computed,
which costs more and answers about the MATRIX: a factorization that
failed by rounding luck on a matrix that is in fact positive definite is
recovered there, and one that is genuinely rank deficient is refused
deterministically rather than according to the platform.

## See also

[`solve_pd`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
[`statmod_marginal`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
