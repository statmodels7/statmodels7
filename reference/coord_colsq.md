# Weighted Sums of Squares of Centered Columns

\\\sum_i w_i (x\_{ik} - m_k)^2\\ for the columns \\k\\ asked for, which
is the curvature a coordinate of a centered block has.

## Usage

``` r
coord_colsq(X, w, k, m)
```

## Arguments

- X:

  The block, dense or `dgCMatrix`.

- w:

  The working weights, length `nrow(X)`.

- k:

  The columns, one-based.

- m:

  The weighted means of every column of the block.

## Value

A numeric vector, one entry per column in `k`.

## Details

It is computed as the sum it is, not as \\\sum_i w_i x\_{ik}^2 - W
m_k^2\\, which is the same number in exact arithmetic and loses as many
digits as the mean exceeds the spread. On a sparse block the sum runs
over the stored entries and adds \\m_k^2\\ times the weight of the rows
the column does not store, which is what a stored zero would contribute,
so the block is never densified.

A column that is constant over the rows is the intercept itself once
centered, and its sum is rounding. Where it falls below \\10^{-10}\\ of
the uncentered sum the uncentered one is returned instead: the kernel
leaves such a coordinate where it is, and the value here only has to be
a usable step for the table.

## See also

[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md),
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
