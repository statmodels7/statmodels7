# Bind a Model's Term Blocks Side by Side

The assembled design, sparse when any term's block is.

## Usage

``` r
bind_blocks(mats, n)
```

## Arguments

- mats:

  A list of blocks.

- n:

  The number of observations, for the empty case.

## Value

A matrix, or a Matrix sparse matrix.

## Details

A grouping indicator is sparse by construction – a row belongs to one
group, so a random effect over \\m\\ of them has density \\1/m\\ – and
modelterms7 builds it that way. Binding it beside a dense block with
[`cbind()`](https://rdrr.io/r/base/cbind.html) does not work: base
dispatch reads the sparse block as a vector and reports that the number
of items to replace is not a multiple of the replacement length, three
frames from anything a caller wrote.

The result is sparse whenever ONE block is, which is the right rule
because sparsity is a property of the assembled matrix rather than of
its pieces: a dense fixed block beside a large indicator leaves a matrix
that is still overwhelmingly zero, and its factorization stays sparse
under a fill-reducing ordering.
