# Zero the Non-Finite Entries of a Penalty's Hessian

A penalty at a hyperparameter far enough out returns non-finite entries,
and every consumer of
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
zeroes them before using the matrix.

## Usage

``` r
zap_nonfinite(S)
```

## Arguments

- S:

  A penalty's Hessian, a square matrix, sparse or dense.

## Value

`S` with its non-finite entries replaced by zero, in the storage it
arrived in.

## Details

Seven places wrote `S[!is.finite(S)] <- 0`. That is correct on a base
matrix and a trap on a sparse one: the logical index is a dense \\p
\times p\\ matrix, so the storage the accumulator exists to keep is
thrown away at the first consumer.

On a sparse matrix only the **stored** values can be non-finite, a
structural zero being finite by construction, so the same answer comes
from the value slot alone in \\O(\mathrm{nnz})\\. The two branches are
written once here in place of seven times.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
