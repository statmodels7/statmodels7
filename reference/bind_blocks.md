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
  there is no block from which to read a row count.

## Value

An `n x p` matrix with `p` the total width of the blocks: a Matrix
object when any input is sparse, a base matrix otherwise. For an empty
`mats` the result is an `n x 0` base matrix, which is the design of an
equation that carries only a structural term.

## Details

A grouping indicator is sparse by construction. A row belongs to one
group, so a random effect over \\m\\ groups has density \\1/m\\, and
modelterms7 builds it that way.

Blocks are bound with [`cbind()`](https://rdrr.io/r/base/cbind.html)
when all of them are base matrices. When any block is sparse, they are
bound with [`Matrix::cbind2()`](https://rdrr.io/r/methods/cbind2.html),
so that the result is a Matrix object whatever the order of the blocks.

The result is sparse when one block is. Sparsity is a property of the
assembled matrix and not of its pieces: a dense fixed block beside a
large indicator leaves a matrix that is still overwhelmingly zero, and
its factorization stays sparse under a fill-reducing ordering.

## See also

[`design_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
for the same test applied to a whole design,
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for the assembly that this serves.
