# The Matrix the Traces Are Taken Against

\\M\\, which is \\K^{-1}\\ when nothing is projected away and the
projected inverse otherwise.

## Usage

``` r
ctx_trace_matrix(ctx, pen, basis)
```

## Arguments

- ctx:

  A context, or `NULL`.

- pen:

  The result of
  [`ctx_penalized`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).

- basis:

  The integrated subspace, or `NULL`.

## Value

A square matrix, or `NULL`.
