# An Orthonormal Basis of a Penalty's Range Space

The directions a penalty constrains, which are the ones
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates over.

## Usage

``` r
penalty_range_basis(pen, k, p, nm)
```

## Arguments

- pen:

  A penalties7 penalty.

- k:

  The number of coefficients in the term's block.

- p:

  The distribution parameter, for the message.

- nm:

  The term's name, for the message.

## Value

A `k` by `r` matrix.
