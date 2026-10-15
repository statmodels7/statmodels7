# When a Matrix Is Factorized Sparsely

Whether a symmetric matrix is large enough and sparse enough that a
sparse Cholesky factorization is used in place of a dense one. The
decision depends on the matrix alone: no term, no model and no caller's
preference enters.

## Usage

``` r
worth_sparse(M, min_dim = 100L, max_density = 0.1)
```

## Arguments

- M:

  A symmetric matrix.

- min_dim:

  The smallest order for which a sparse factorization is used.

- max_density:

  The largest fraction of nonzero entries for which a sparse
  factorization is used.

## Value

A single logical.

## Details

A matrix is factorized sparsely when it is stored as an S4 sparse
matrix, its order is at least `min_dim` and its fraction of nonzero
entries is at most `max_density`. Below about a hundred coefficients the
sparse route is slower, because its fixed cost is the coercion and the
S4 dispatch around it, which does not shrink with the matrix. The size
condition prevents this on the fully dense penalized information of a
single smooth.

**Both quantities are read from the matrix, and the first is its
storage.** A matrix held as a base matrix is never factorized sparsely,
whatever its zeros. This is not a test of the term that produced the
matrix.
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
accumulates into the design's own kind, so the penalized matrix is
stored sparsely exactly when the design is. modelterms7 settles the
storage of a block at build when `sparse = NULL`, the default, and
builds it sparse where the product of the number of rows and the number
of indicator columns exceeds \\10^5\\; `sparse = TRUE` and
`sparse = FALSE` force the sparse and the dense storage. Where the
storage is dense the factorization is a small share of the time of a
fit, which is spent in the \\O(np^2)\\ products against the design, and
coercing a dense \\p \times p\\ matrix here to save a share of that size
would cost more than it returns.

## See also

[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
