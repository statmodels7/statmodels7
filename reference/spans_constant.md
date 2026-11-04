# Whether a Design Spans the Constant

Returns `TRUE` when regressing the constant on the columns of a design
leaves no residual.

## Usage

``` r
spans_constant(X, n)
```

## Arguments

- X:

  A matrix, possibly with no column.

- n:

  The number of rows.

## Value

A single logical.

## See also

[`unpenalized_design()`](https://statmodels7.github.io/statmodels7/reference/unpenalized_design.md)
