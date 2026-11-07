# The Coordinates of a Pinned Matrix That Are Not Pinned

Returns the coordinates of a matrix that
[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md)
did not set to the identity, or `NULL` where nothing is pinned or the
pinned rows are not decoupled from the rest.

## Usage

``` r
pinned_free(M)
```

## Arguments

- M:

  A matrix, dense or sparse, possibly carrying the attribute `held_at`.

## Value

An integer vector, possibly empty, or `NULL`.

## Details

A pinned coordinate is a row and a column of the identity, so the matrix
is the block of the other coordinates beside an identity block, and its
smallest eigenvalue is at most 1. The condition test of
[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
compares the smallest eigenvalue with \\\sqrt{\epsilon}\\ times the norm
of the matrix, so on the whole matrix it rejected every Cholesky factor
once the norm passed \\1/\sqrt{\epsilon}\\, whatever the conditioning of
the block. The test is therefore read on the block of the coordinates
returned here. The rows named by the attribute `held_at` are checked to
be rows of the identity, so a matrix that was changed after it was
pinned keeps the test on the whole matrix.

## See also

[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md),
[`pd_logdet_dense()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet_dense.md),
[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md)
