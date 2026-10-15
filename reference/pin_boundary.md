# Pin the Coordinates a Boundary Has Frozen

Replaces the row and column of every coordinate whose curvature is not
finite, and of every coordinate named in `hold`, by the identity's,
leaving the rest of the matrix alone.

## Usage

``` r
pin_boundary(K, hold = integer(0))
```

## Arguments

- K:

  The penalized information.

- hold:

  Integer positions to pin in addition to the coordinates with a
  non-finite diagonal (the coefficients that the marginal criterion
  holds, and the coordinates that a kinked penalty covers). None by
  default.

## Value

`K` with the pinned coordinates replaced by the identity's rows and
columns, carrying the attributes `held` (their number), `held_at` (their
positions), `boundary_at` (the positions with a non-finite diagonal) and
`hold_at` (the positions pinned through `hold` alone). Where there are
no pinned coordinates the matrix is otherwise unchanged and the
positions are empty.

## Details

A parameter that has reached the clamp that its link keeps it strictly
inside leaves the family's curvature there at `NaN`, and one such entry
is enough to deny the whole matrix a Cholesky factor. The coordinate is
not one that the criterion integrates over, being held at a boundary
where its own score is exactly zero, so it contributes nothing to the
Laplace approximation. A unit diagonal contributes exactly that:
\\\log\lvert K\rvert\\ gains \\\log 1 = 0\\, and \\K^{-1}g\\ returns
that coordinate's own score, which is zero.

The dimension of `K` is preserved, because its consumers index `K` and
its inverse by position. Dropping the coordinate would change both
dimensions.

## See also

[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md),
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
