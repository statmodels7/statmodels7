# Zero the Non-Finite Entries of a Penalty's Hessian

A penalty at a hyperparameter far enough out returns non-finite entries,
and every consumer of
[`statmod_penalty_at`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
zeroes them before using the matrix.

## Usage

``` r
zap_nonfinite(S)
```

## Arguments

- S:

  A penalty's Hessian, sparse or dense.

## Value

`S` with its non-finite entries replaced by zero.

## Details

Seven places wrote `S[!is.finite(S)] <- 0`, which is correct on a base
matrix and a trap on a sparse one: the logical index is a DENSE \\p
\times p\\ matrix, so the storage the accumulator exists to keep is
thrown away at the first consumer. On a sparse matrix only the STORED
values can be non-finite – a structural zero is finite by construction –
so the same answer comes from the value slot alone, in
\\O(\mathrm{nnz})\\.

Written once here rather than seven times, which is the discipline this
package records for a shape of mistake rather than a line of it.

## See also

[`statmod_penalty_at`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
