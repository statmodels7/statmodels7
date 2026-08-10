# The Square-Root Design

Returns \\R\\ with \\R'R = Z'\Omega Z\\, the matrix a scoring step
decomposes instead of the information itself.

## Usage

``` r
sqrt_design(design, L)
```

## Arguments

- design:

  The design, as
  [`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- L:

  The Cholesky factors, from
  [`chol_blocks`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md).

## Value

An \\nK \times p\\ matrix, or `NULL` when a block was not positive
definite.

## Details

Row block \\a\\ of \\R\\ carries, in the columns of parameter \\b\\, the
entry \\L_i\[b, a\]\\ times that parameter's design row; the factor
being lower triangular, the blocks with \\b \< a\\ are zero and are not
formed.
