# Solve a Scoring Step From the Square-Root Design

Decomposes the augmented matrix \\A = \[R;\\ C\]\\, whose cross-product
\\A'A = R'R + C'C\\ is the penalized information, and returns the
increment solving

\$\$(R'R + C'C)\\\delta = u\$\$

The augmented form is the point: neither \\R'R\\ nor the penalized
information is ever formed, so a design whose condition number is
\\\kappa\\ is decomposed at \\\kappa\\. Forming the normal equations
would square it.

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

  `"qr"` or `"svd"`. `"svd"` has no sparse counterpart, so a sparse pair
  asked for it is densified first.

- threads:

  How many threads the triangular factor may use, a plain integer. `1L`
  takes [`qr()`](https://rdrr.io/r/base/qr.html) unconditionally.

## Value

A list of two:

- `delta`:

  the increment, an unnamed numeric vector of length `p`. Zero in any
  coordinate the pivoted route dropped.

- `rank`:

  the rank used, an integer. Always `p` on the threaded and sparse
  routes, which decline on a deficient matrix instead of reporting a
  smaller rank.

## The four routes

Tried in this order.

1.  **Sparse QR**, when either factor is a Matrix object. Handed to
    [`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md),
    which declines on a rank-deficient matrix.

2.  **SVD**, when `how = "svd"`. The singular values below \\\max(\dim
    A)\\\epsilon\\\sigma\_{\max}\\ are dropped and the increment is
    built from the retained right singular vectors, so a deficient
    system gets the minimum-norm answer.

3.  **The threaded triangular factor**, when `how = "qr"`, `threads > 1`
    and the work \\n p^2\\ is at least `5e7`. Only the triangular factor
    is ever read, so a kernel that produces it and never accumulates
    \\Q\\ does the whole job.

4.  **[`qr()`](https://rdrr.io/r/base/qr.html)**, the pivoted LINPACK
    route, which reports a rank and can drop columns. A sparse solve
    that declined falls through to here on densified factors.

## Why the rank test is equilibrated

The threaded route is taken only where the matrix is comfortably of full
rank, and the test reads the diagonal of the triangular factor **divided
by each column's norm**. Since \\A'A = R'R\\, scaling \\A\\'s columns
scales that diagonal by the same factors, so the ratio is the diagonal
of a decomposition with unit column norms.

Reading the raw diagonal instead reports a matrix as near-singular
whenever its columns differ in size, and a large smoothing parameter
does exactly that to the penalty rows of its own block beside an
unpenalized one. Per-column scaling forgives separation from any source;
an exact collinearity stays exactly singular.

The tolerance is `1e-7` on that ratio, the same one `dqrdc2` uses, so
anything the pivoted route would call deficient still goes there and is
reported with its rank.

Measured against [`qr()`](https://rdrr.io/r/base/qr.html) at \\n =
40000\\ with the column scales spread over \\10^8\\: 2.6x at \\p = 51\\
and 3.8x at \\p = 145\\ to 600 on eight threads, with the increment
agreeing to \\1.6 \times 10^{-16}\\.

## See also

[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
and
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
for the two factors,
[`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md)
for the sparse route.
