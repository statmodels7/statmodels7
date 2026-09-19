# A Symmetric Matrix From Its Upper Blocks

Binds a symmetric matrix out of the blocks on and above its block
diagonal, the ones below being the transposes of those above.

## Usage

``` r
assemble_blocks(blocks, sparse)
```

## Arguments

- blocks:

  A list of lists, `blocks[[a]][[b]]` the block of rows `a` and columns
  `b` for `b >= a`. Entries below the diagonal are not read.

- sparse:

  Whether the result is a `dgCMatrix`, as
  [`design_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  decides for the accumulator it replaces.

## Value

The assembled matrix, without dimnames: a `dgCMatrix` when `sparse` is
`TRUE`, otherwise a base matrix.

## Details

It replaces writing each block into a zero accumulator. With a sparse
design that accumulator is a `dgCMatrix`, and a sub-assignment into it
rebuilds the compressed columns every time: measured on a lasso over 200
columns beside a random intercept over 500 groups, the three assignments
of one information cost 0.07 s against 0.19 s for the products
themselves. Binding the blocks once does no arithmetic on them, so the
entries are the ones the products gave.

## See also

[`zero_information()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md),
the accumulator this replaces where every block is written exactly once.
