# The Column Rank of One Design Block

The number of linearly independent columns, by the pivoted decomposition
the fit's own solve uses, so that the two cannot disagree about which
columns a model identifies.

## Usage

``` r
block_column_rank(X)
```

## Arguments

- X:

  A design block, dense or sparse.

## Value

A single number, or `NA_integer_` where the decomposition fails.

## Details

A dense block goes through [`qr()`](https://rdrr.io/r/base/qr.html), at
`dqrdc2`'s own tolerance, which is what
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
reads. A sparse one goes through
[`Matrix::qr()`](https://rdrr.io/pkg/Matrix/man/qr-methods.html) and its
rank is counted on the JACOBI-EQUILIBRATED diagonal of the triangular
factor, which is the correction
[`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md)
already carries: since \\R'R = A'A\\, scaling the columns by their norms
scales that diagonal by the same factors, so a block whose columns
differ in size is not read as deficient while an exact dependence stays
exactly singular.

## See also

[`design_count_exact()`](https://statmodels7.github.io/statmodels7/reference/design_count_exact.md),
[`active_rank()`](https://statmodels7.github.io/statmodels7/reference/active_rank.md).
