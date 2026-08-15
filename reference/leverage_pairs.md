# The Leverage Diagonal Over the Nonzeros of Two Rows

\\G_i = \sum\_{j\in J_i}\sum\_{k\in K_i}
X_a\[i,j\]X_b\[i,k\]M\_{ab}\[j,k\]\\, the same quantity
[`block_leverage`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
otherwise reads off a dense \\n\times p_b\\ product.

## Usage

``` r
leverage_pairs(ta, tb, Mab, n)
```

## Arguments

- ta, tb:

  The two designs' nonzeros, from
  [`row_nonzeros`](https://statmodels7.github.io/statmodels7/reference/row_nonzeros.md).

- Mab:

  The block of \\M\\.

- n:

  The number of observations.

## Value

A numeric vector as long as the sample.

## Details

Where a design is built from grouping indicators a row has one nonzero
per block, so the quadratic form is over a handful of entries and the
dense product computes \\p_a p_b\\ of them per row to keep one. The
pairs are expanded once and the whole sum is vectorized.

**It is taken only where it wins**, and the threshold is measured rather
than assumed. At a combined density of 3.6e-05 (a random intercept over
500 levels) it is 14.2 times the dense route; at 0.18 (three smooths and
a random effect) it is 50 times SLOWER, R's per-element indexing being
far dearer than a BLAS flop, and at a dense block 50 times slower again.
Interpolating the two measurements puts the crossover at a combined
density near 1.1e-03, which is the gate: below it the route is taken, at
it the two cost the same, and above it the dense product stands.

## See also

[`block_leverage`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
