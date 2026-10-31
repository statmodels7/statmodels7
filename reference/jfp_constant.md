# Whether the Constant Direction Enters a Jeffreys Block

Applies the rule of the `intercept` argument of `jfp()` to the columns
of a formula level and to the columns of the block.

## Usage

``` r
jfp_constant(intercept, X_level, X_A, tol = 1e-07)
```

## Arguments

- intercept:

  `NULL`, `TRUE` or `FALSE`.

- X_level:

  The columns of the formula level outside the block, or `NULL`.

- X_A:

  The block's columns.

- tol:

  The tolerance of the span test.

## Value

A single logical.

## Details

With `intercept = FALSE` the direction is never included. Where the
column of ones lies in the span of the block's own columns the direction
is not included either, whatever `intercept` says, because the
determinant would be zero and the penalty is already invariant to
centering a column. Otherwise `TRUE` includes it, and `NULL` includes it
when the column of ones lies in the span of `X_level`, that is when the
model is invariant to centering a column of the block.

The span is tested by
[`ones_in_span()`](https://statmodels7.github.io/statmodels7/reference/ones_in_span.md).

## See also

[`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
