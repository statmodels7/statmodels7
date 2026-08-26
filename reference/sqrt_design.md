# The Square-Root Design

Builds \\R\\ with \\R'R = Z'\Omega Z\\, the square root of the assembled
information. A scoring step decomposes this instead of the information
itself, so \\Z'\Omega Z\\ is never formed and its condition number is
never squared.

## Usage

``` r
sqrt_design(design, L)
```

## Arguments

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it, holding one block of columns per distribution parameter.

- L:

  The per-observation Cholesky factors, as
  [`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
  returns them, or `NULL`.

## Value

An \\nK \times p\\ matrix, dense or a Matrix sparse matrix according to
the design's own storage, where \\p\\ is the total number of
coefficients across the equations. `NULL` when `L` is `NULL`, so the
refusal propagates from the factorization to the solve without a test at
each step.

## Details

Row block \\a\\ of \\R\\ carries, in the columns belonging to parameter
\\b\\, the entry \\L_i\[b, a\]\\ times that parameter's design row. The
factor is lower triangular, so the blocks with \\b \< a\\ are zero and
are not formed at all.

The result is sparse when any equation's design block is, and the solve
then goes to
[`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md).

## See also

[`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
for `L`,
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
for the decomposition this feeds,
[`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
for the other half of the augmented matrix.
