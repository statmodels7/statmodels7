# The Augmented Solve With a Product of Low Rank

Solves \\(R^\top R + C^\top C + VV^\top)\delta = u\\ from one
decomposition of \\\[R; C\]\\, by Woodbury's identity.

## Usage

``` r
augmented_lowrank(R, C, V, u, how, threads = 1L, tol = NULL)
```

## Arguments

- R, C:

  The square-root design and the penalty's factor.

- V:

  The low-rank factor, or `NULL` for the plain augmented solve.

- u:

  The right-hand side.

- how, threads, tol:

  As in
  [`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md).

## Value

What
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
returns.

## Details

[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
solves the right-hand side and the columns of \\V\\ with one
factorization, and the correction is a system of the size of \\V\\'s
columns. Appending \\V^\top\\ as rows of the augmented matrix would give
the same solution, but each such row is dense, and a sparse QR
decomposition fills in completely beside it: measured on a random
intercept over 250 levels, the Jeffreys fit took 13.9 s against 3.1 s
without the block, almost all of it in that decomposition. A coordinate
that the pivot leaves out is read as zero in every column.

## See also

[`lowrank_solve()`](https://statmodels7.github.io/statmodels7/reference/lowrank_solve.md),
the same identity on the assembled route.
