# The Leverage Diagonal Over the Nonzeros of Two Rows

Computes \\G_i = \sum\_{j\in J_i}\sum\_{k\in K_i}
X_a\[i,j\]X_b\[i,k\]M\_{ab}\[j,k\]\\, the same quantity that
[`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
otherwise reads off a dense \\n\times p_b\\ product.

## Usage

``` r
leverage_pairs(ta, tb, Mab, n, threads = 1L)
```

## Arguments

- ta, tb:

  The two designs' nonzeros, from
  [`row_nonzeros()`](https://statmodels7.github.io/statmodels7/reference/row_nonzeros.md).

- Mab:

  The block of \\M\\.

- n:

  The number of observations.

- threads:

  How many threads the kernel may use.

## Value

A numeric vector as long as the sample.

## Details

Where a design is built from grouping indicators a row has one nonzero
per block, so the quadratic form is over a handful of entries, while the
dense product computes \\p_a p_b\\ of them per row to keep one. The
compiled kernel `leverage_pairs_cpp()` reads the nonzeros of the two
designs row by row and sums them without expanding the pairs, and each
row is written by one thread.
[`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
takes this route only where both blocks are sparse, the product of their
densities is below \\10^{-3}\\ and the number of pairs is below \\2
\times 10^{7}\\.

## See also

[`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
