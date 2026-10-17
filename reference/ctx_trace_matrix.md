# The Matrix Against Which the Traces Are Taken

Returns \\M\\, which is the projected inverse when a subspace is
projected away and \\K^{-1}\\ otherwise.

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

  The information with which `pen` was built, which is the cache's key:
  the projection is of that matrix, so one entry cannot serve both
  information types.

## Value

A square matrix in which the rows and columns of the pinned coordinates
are zero, or `NULL`.
