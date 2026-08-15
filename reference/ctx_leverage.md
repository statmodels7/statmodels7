# The Per-Observation Diagonals at the Context's Point

[`block_leverage`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md),
computed once and read by the gradient's contraction and by every pair
of the Hessian.

## Usage

``` r
ctx_leverage(ctx, design, M, params, npar, offs)
```

## Arguments

- ctx:

  A context, or `NULL`.

- design:

  The design.

- M:

  The matrix the traces are taken against.

- params, npar, offs:

  The block bookkeeping.

## Value

A list of lists of vectors.
