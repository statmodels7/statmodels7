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

`S`, in which the non-finite entries are replaced by zero, and which
keeps the storage in which it arrived.

## Details

The assignment `S[!is.finite(S)] <- 0` is correct on a base matrix. On a
sparse matrix the logical index is a dense \\p \times p\\ matrix, so the
sparse storage that the accumulator is meant to preserve would be lost.

On a sparse matrix only the stored values can be non-finite, because a
structural zero is finite by construction, so the same result comes from
the value slot alone in \\O(\mathrm{nnz})\\ operations.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
