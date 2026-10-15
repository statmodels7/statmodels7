# Factorize a Penalized Information Once

The Cholesky factor of a matrix that a Laplace approximation requires to
be positive definite, together with its log-determinant, in the storage
of the matrix itself.

## Usage

``` r
pd_factor(M, scale = NULL)
```

## Arguments

- M:

  A symmetric matrix, sparse or dense.

- scale:

  A reference magnitude, as
  [`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
  takes.

## Value

A list with `logdet`, `ok`, `factor` and `sparse`. The factor is `NULL`
where the answer came from the eigendecomposition.

## Details

This function is the single place where the penalized matrix is
factorized. The criterion needs its log-determinant, the gradient needs
the movement of the mode and the Hessian needs both and the inverse, so
one factorization serves all three.

**Definiteness.** Whether the matrix is accepted never depends on
whether a factorization raised an error: where the cheap test is
inconclusive, the eigendecomposition answers about the matrix. The
sparse route carries its own condition estimate
([`sparse_lmin()`](https://statmodels7.github.io/statmodels7/reference/sparse_lmin.md))
in place of the dense one, and falls back to the dense route where that
estimate cannot be formed, so the matrix is judged by the same criterion
on either storage.

## See also

[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md),
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
