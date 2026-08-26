# Solve a Scoring Step From a Sparse Square-Root Design

The same increment
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
returns, taken through a sparse QR of \\\[R;\\ C\]\\.

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

  The decomposition asked for. Only `"qr"` is served; `"svd"` has no
  sparse counterpart and declines.

## Value

A list with `delta` (numeric, length `p`) and `rank` (integer, equal to
`p`), or `NULL` when the route declines: `how` is not `"qr"`, the
factorization fails, or the equilibrated diagonal says the matrix is
deficient.

## Why a QR

The augmented form exists so that \\X'X\\ is never formed and the
conditioning is never squared. A sparse QR is a QR and keeps that
property exactly. A sparse Cholesky of the normal equations would be
faster still and would give it up.

Measured against the dense QR on the same augmented design of a
random-intercept model: 695 times faster at 100 groups and 75475 times
at 1000, where the dense factorization costs 9.06 s against 0.00012 s.

## The factor is not back-permuted

A sparse QR reorders the columns to reduce fill, so \\AP = QR\\ for the
permutation \\P\\ the decomposition chose, and \\(A'A)^{-1} =
P(R'R)^{-1}P'\\. The increment is two triangular solves between a
permutation and its inverse, which is the bookkeeping the dense route
already does with [`qr()`](https://rdrr.io/r/base/qr.html)'s pivot.

`qrR(backPermute = TRUE)` looks like the simpler choice and is a trap.
The back-permuted factor is no longer triangular, so its diagonal says
nothing about the rank, and a solve against it is a general solve
instead of two triangular ones.

## The rank test

Read off the diagonal of the triangular factor, since a sparse QR does
not report a rank. The test is on the **Jacobi-equilibrated** diagonal,
each entry divided by its column's norm, and declines when the smallest
ratio falls to `ncol(A) * .Machine$double.eps` of the largest.

Without the equilibration the route is barely usable. Reading the raw
diagonal rejected 87 of 127 solves on matrices the dense route finds at
full rank, the ratio there running down to \\7 \times 10^{-30}\\ while
the equilibrated one stayed between 0.445 and 1. Those matrices are not
near-singular; their columns differ in size, as a large smoothing
parameter makes them. With the raw test a random-effect fit fell through
to a dense QR and cost 34.1 s where it now costs 3.9.

The column norms come from `colSums(A^2)`. Writing `colSums(A * A)`
instead makes a binary operation between two sparse matrices, which
intersects their index sets and is 29 to 35 times slower; the norms were
70 to 74 per cent of the whole solve until that line changed.

Where the matrix really is rank deficient there is no unique increment,
and this route declines instead of choosing one. The caller falls back
to the dense route, which drops columns and says how many it kept.

## See also

[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md),
the caller and the dense fallback.
