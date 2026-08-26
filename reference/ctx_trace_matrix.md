# The Matrix the Traces Are Taken Against

\\M\\, which is \\K^{-1}\\ when nothing is projected away and the
projected inverse otherwise.

## Usage

``` r
ctx_trace_matrix(ctx, pen, basis, expected = FALSE)
```

## Arguments

- ctx:

  A context, or `NULL`.

- pen:

  The result of
  [`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).

- basis:

  The integrated subspace, or `NULL`.

- expected:

  Which information `pen` was built with, which is the cache's key: the
  projection is of that matrix, so one entry cannot serve both. Nothing
  reaches it today, a search holding one
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  throughout. The slot is keyed all the same, the twin defect in
  [`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
  having been unreachable in exactly the same way until it was not.

## Value

A square matrix, or `NULL`.
