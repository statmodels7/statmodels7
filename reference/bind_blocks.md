# Bind a Model's Term Blocks Side by Side

Binds one equation's term blocks side by side into its design, keeping
the result sparse whenever any block is.

## Usage

``` r
bind_blocks(mats, n)
```

## Arguments

- mats:

  A list of design blocks, all with the same number of rows, each dense
  or sparse. May be empty.

- n:

  The number of observations, used only when `mats` is empty, where
  there is no block to read a row count from.

## Value

An `n x p` matrix with `p` the total width of the blocks: a Matrix
object when any input is sparse, a base matrix otherwise. An `n x 0`
base matrix for an empty `mats`, as an equation carrying only a
structural term gives.

## Details

A grouping indicator is sparse by construction. A row belongs to one
group, so a random effect over \\m\\ groups has density \\1/m\\, and
modelterms7 builds it that way.

Binding such a block beside a dense one with
[`cbind()`](https://rdrr.io/r/base/cbind.html) fails: base dispatch
reads the sparse block as a vector and reports that the number of items
to replace is not a multiple of the replacement length, three frames
from anything a caller wrote.
[`Matrix::cbind2()`](https://rdrr.io/r/methods/cbind2.html) is what
handles the mixed pair.

The result is sparse when **one** block is. Sparsity is a property of
the assembled matrix, never of its pieces: a dense fixed block beside a
large indicator leaves a matrix that is still overwhelmingly zero, and
its factorization stays sparse under a fill-reducing ordering.

## See also

[`design_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
for the same question asked of a whole design,
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for the assembly this serves.
