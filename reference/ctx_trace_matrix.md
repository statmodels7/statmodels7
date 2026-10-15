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
  projection is of that matrix, so one entry cannot serve both
  information types.

## Value

A square matrix, with the rows and columns of the pinned coordinates set
to zero, or `NULL`.
