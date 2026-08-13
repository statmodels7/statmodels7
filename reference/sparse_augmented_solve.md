# Solve a Scoring Step From a Sparse Square-Root Design

The same increment
[`augmented_solve`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
returns, taken through a sparse QR of \\\[R;\\ C\]\\.

## Usage

``` r
sparse_augmented_solve(R, C, u, how)
```

## Arguments

- R:

  The square-root design.

- C:

  The penalty's factor.

- u:

  The right-hand side.

- how:

  The decomposition asked for; `"svd"` has no sparse counterpart and
  declines.

## Value

A list with `delta` and `rank`, or `NULL`.

## Details

The augmented route exists so that \\X'X\\ is never formed and the
conditioning is never squared, and a sparse QR is a QR: it keeps that
property exactly, which a sparse Cholesky of the normal equations would
not. Measured against the dense QR on the same augmented design of a
random-intercept model, it is 695 times faster at 100 groups and 75475
at 1000, the dense factorization there costing 9.06 s against 0.00012.

The factor is taken WITHOUT back-permuting. A sparse QR reorders the
columns to reduce fill, so \\AP = QR\\ for the permutation \\P\\ the
decomposition chose, and \\(A'A)^{-1} = P(R'R)^{-1}P'\\: the increment
is two triangular solves between a permutation and its inverse, which is
the same bookkeeping the dense route does with `qr`'s pivot.
`backPermute = TRUE` looks simpler and is a trap – it returns a factor
that is no longer triangular, so its diagonal says nothing about the
rank and a solve against it is a general one rather than two triangular
ones.

The rank is read off the diagonal of the triangular factor rather than
reported by the decomposition, which for a sparse QR does not give one.
Where the matrix is rank deficient there is no unique increment to
return and this route declines, leaving the caller its dense fallback,
which can drop columns and say how many it kept.

## See also

[`augmented_solve`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
