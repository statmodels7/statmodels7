# The Derivatives of a Refreshed Jeffreys Block's Rows

`jfp_vk()` returns the derivative of the block's rows in one of the
term's coefficients, an observations by \\\tilde p\\ matrix, read from
the state's `D`.

## Usage

``` r
jfp_vk(bs, l)
```

## Arguments

- bs:

  A block's entry of the state, of kind `"refresh"`.

- l:

  A position among the term's coefficients.

## Value

A numeric matrix.
