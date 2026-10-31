# Whether the Column of Ones Lies in the Span of a Matrix

Compares the rank of the matrix with the rank of the matrix with a
column of ones appended, both from the pivoted QR decomposition with the
tolerance with which the toolkit detects aliased coefficients.

## Usage

``` r
ones_in_span(X, tol = 1e-07)
```

## Arguments

- X:

  A matrix, a Matrix object, or `NULL`.

- tol:

  The tolerance of [`qr()`](https://rdrr.io/r/base/qr.html).

## Value

A single logical.
