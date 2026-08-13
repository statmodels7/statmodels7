# The Penalized Block, in the Storage It Arrived In

Slices a penalized term's columns out of its equation's design without
densifying a sparse one, and normalizes a Matrix to the
compressed-column class the kernel reads.

## Usage

``` r
coord_block(X, cols)
```

## Arguments

- X:

  The equation's design.

- cols:

  The term's column positions.

## Value

A numeric matrix or a `dgCMatrix`.

## Details

A dense slice of a base matrix is returned as it is. Any Matrix is
carried to `dgCMatrix`: the general compressed-column form is the one
whose slots the kernel walks, and a symmetric or triangular compression
would describe the same entries differently. A dense Matrix class is
materialized as a base matrix instead, there being nothing to save.

## See also

[`coord_call`](https://statmodels7.github.io/statmodels7/reference/coord_call.md),
[`coord_fit`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
