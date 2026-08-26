# Factorize a Penalized Information Once

The Cholesky factor of a matrix a Laplace approximation needs positive
definite, together with its log-determinant, in whichever storage the
matrix itself calls for.

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

This is the one place the penalized matrix is factorized. The criterion
wants its log-determinant, the gradient wants the mode's movement and
the Hessian wants both plus the inverse; before this existed the
criterion and
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
each factorized the same matrix at the same point, which at \\p = 503\\
was 12.4 ms spent twice.

**The verdict is unchanged and so is its property.** Whether the matrix
is accepted never turns on whether a factorization raised: where the
cheap test is inconclusive the eigendecomposition answers about the
matrix. The sparse route carries its own condition estimate
([`sparse_lmin()`](https://statmodels7.github.io/statmodels7/reference/sparse_lmin.md))
rather than the dense one, and falls back to the dense route where that
estimate cannot be formed, so a refusal is reached by the same reasoning
on either storage.

## See also

[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md),
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
