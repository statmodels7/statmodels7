# The Penalized Matrix With a Jeffreys Curvature, Factored in Two Parts

Returns what
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
returns, for a model carrying a Jeffreys term, from a factorization of
the sparse part of \\K\\ and Woodbury's identity for its part of low
rank, or `NULL` where the model carries no Jeffreys term or the sparse
part has no positive definite factor, and
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
then assembles \\K\\ whole.

## Usage

``` r
penalized_lowrank(ctx, spec, design, coef, hyper, H, hold)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  As in
  [`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).

- H:

  The information.

- hold:

  The coordinates to pin, beside the boundary ones.

## Value

A list with `K`, `inv` and `logdet`, or `NULL`.

## Details

With \\K = K_0 + VV^\top\\, \$\$\log\|K\| = \log\|K_0\| + \log\|I +
V^\top K_0^{-1}V\|,\qquad K^{-1} = K_0^{-1} - K_0^{-1}V(I + V^\top
K_0^{-1}V)^{-1}V^\top K_0^{-1}.\$\$ A coordinate held in \\K_0\\ has a
zero row in \\V\\, which is the row that pinning the assembled matrix
would leave. The matrix returned as `K` is the assembled sum, for its
readers that multiply by it.

## See also

[`lowrank_solve()`](https://statmodels7.github.io/statmodels7/reference/lowrank_solve.md),
the same identity in the inner fit.
