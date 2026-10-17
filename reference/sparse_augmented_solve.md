# Solve a Scoring Step From a Sparse Square-Root Design

Returns the same increment that
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
returns, computed through a sparse QR of \\\[R;\\ C\]\\.

## Usage

``` r
sparse_augmented_solve(R, C, u, how)
```

## Arguments

- R:

  The square-root design, \\nK \times p\\, sparse or dense.

- C:

  The penalty's factor, with `p` columns.

- u:

  The right-hand side, a numeric vector of length `p`.

- how:

  The decomposition requested. Only `"qr"` is served; `"svd"` has no
  sparse counterpart, so the function returns `NULL` for it.

## Value

A list with `delta` (numeric, length `p`) and `rank` (integer, equal to
`p`), or `NULL` when the route yields no solution: `how` is not `"qr"`,
the stacked matrix has fewer rows than columns, a column has zero norm,
the factorization or the triangular solves fail, or the equilibrated
diagonal shows that the matrix is deficient. Unlike
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md),
the list has no `dropped` element.

## Why a QR

The augmented form exists so that \\X'X\\ is never formed and the
conditioning is never squared. A sparse QR is a QR and keeps that
property. A sparse Cholesky factorization of the normal equations would
be faster but would lose that property.

## The factor is not back-permuted

A sparse QR reorders the columns to reduce fill, so \\AP = QR\\ for the
permutation \\P\\ that the decomposition selected, and \\(A'A)^{-1} =
P(R'R)^{-1}P'\\. The increment is two triangular solves between a
permutation and its inverse, which is the bookkeeping that the dense
route does with the pivot of [`qr()`](https://rdrr.io/r/base/qr.html).

The function calls `qrR(backPermute = FALSE)`. The back-permuted factor
is not triangular, so its diagonal does not give the rank, and a solve
against it is a general solve instead of two triangular ones.

## The rank test

The rank is read off the diagonal of the triangular factor, since a
sparse QR does not report a rank. The test is on the Jacobi-equilibrated
diagonal, in which each entry is divided by its column's norm, and the
route returns `NULL` when the smallest ratio is at most
`ncol(A) * .Machine$double.eps` times the largest.

Without the equilibration the raw diagonal would reject matrices that
the dense route finds at full rank, because their columns differ in size
(a large smoothing parameter makes them so), and the solve would fall
through to the dense QR.

The column norms come from `colSums(A^2)`, which acts on the `x` slot of
the sparse matrix; `colSums(A * A)` would form a binary operation
between two sparse matrices, which intersects their index sets and is
slower.

Where the matrix is rank deficient there is no unique increment, and
this route returns `NULL` instead of choosing one. The caller falls back
to the dense route, which drops columns and reports the rank and the
dropped columns.

## See also

[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md),
the caller and the dense fallback.
