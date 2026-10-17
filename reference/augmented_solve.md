# Solve a Scoring Step From the Square-Root Design

Decomposes the augmented matrix \\A = \[R;\\ C\]\\, whose cross-product
\\A'A = R'R + C'C\\ is the penalized information, and returns the
increment solving

\$\$(R'R + C'C)\\\delta = u\$\$

The augmented form never forms \\R'R\\ or the penalized information, so
a design with condition number \\\kappa\\ is decomposed at \\\kappa\\
and not at \\\kappa^2\\, which the normal equations would give.

## Usage

``` r
augmented_solve(R, C, u, how, threads = 1L)
```

## Arguments

- R:

  The square-root design, as
  [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  returns it, \\nK \times p\\.

- C:

  The penalty's factor, as
  [`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
  returns it, with `p` columns and at most `p` rows. May have zero rows
  for an unpenalized model.

- u:

  The right-hand side, a numeric vector of length `p`.

- how:

  `"qr"` or `"svd"`. `"svd"` has no sparse counterpart, so with
  `how = "svd"` a sparse pair is densified first.

- threads:

  How many threads the triangular factor may use, a plain integer. With
  `1L` the threaded route is not used.

## Value

A list with the elements `delta`, `rank` and `dropped`:

- `delta`:

  the increment, an unnamed numeric vector of length `p`. It is zero in
  every coordinate that the pivoted route dropped.

- `rank`:

  the rank used, an integer. It is always `p` on the threaded and sparse
  routes, which hand a deficient matrix to the dense route instead of
  reporting a smaller rank, and it is the number of retained singular
  values on the `"svd"` route.

- `dropped`:

  the integer positions of the columns that the pivoted
  [`qr()`](https://rdrr.io/r/base/qr.html) route left out (the aliased
  coordinates); `integer(0)` on the sparse, threaded and `"svd"` routes.

## The four routes

The routes are tried in this order.

1.  Sparse QR, when either factor is a Matrix object. The pair is passed
    to
    [`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md),
    which returns `NULL` on a rank-deficient matrix.

2.  SVD, when `how = "svd"`. The singular values below \\\max(\dim
    A)\\\epsilon\\\sigma\_{\max}\\ are dropped and the increment is
    built from the retained right singular vectors, so the solution of a
    deficient system is the minimum-norm one.

3.  The threaded triangular factor, when `how = "qr"`, `threads > 1` and
    \\m p^2\\ is at least `5e7`, where \\m\\ is the number of rows of
    \\A\\. Only the triangular factor is read, so a kernel that produces
    it without accumulating \\Q\\ does the whole job.

4.  [`qr()`](https://rdrr.io/r/base/qr.html), the pivoted LINPACK route,
    which returns a rank and can drop columns. A sparse solve that
    returned `NULL` continues here on densified factors.

## Why the rank test is equilibrated

The threaded route is taken only where the matrix has full rank by a
comfortable margin. The test reads the diagonal of the triangular factor
divided by each column's norm. Writing \\T\\ for the triangular factor
of \\A\\, so that \\A'A = T'T\\, scaling the columns of \\A\\ scales the
diagonal of \\T\\ by the same factors, so the ratio is the diagonal of a
decomposition with unit column norms.

The raw diagonal would report a matrix as near-singular whenever its
columns differ in size, and a large smoothing parameter produces that
difference between the penalty rows of its own block and an unpenalized
block. Scaling each column to unit norm removes differences in column
size from any source and leaves an exact collinearity exactly singular.

The route is taken when the smallest equilibrated entry exceeds `1e-7`
times the largest, the tolerance that `dqrdc2` uses, so a matrix that
the pivoted route would call deficient is solved there and its rank is
returned.

## See also

[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
and
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
for the two factors,
[`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md)
for the sparse route.
