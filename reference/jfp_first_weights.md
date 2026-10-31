# The Weights of the First Derivative of a Jeffreys Block

Returns \\s\_{c,t} =
-\sum\_{q,r}h\_{qr,t}\\\partial\_{\eta_c}E\_{qr,t}\\ for every equation
\\c\\, the per-observation weights that
[`jfp_gradient()`](https://statmodels7.github.io/statmodels7/reference/jfp_gradient.md)
contracts against each equation's design and that
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
contracts against the second derivative of a refreshed block.

## Usage

``` r
jfp_first_weights(st, bs, dE)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- bs:

  A block's entry of the state.

- dE:

  The derivative of the expected information, from
  [`jfp_dexpected()`](https://statmodels7.github.io/statmodels7/reference/jfp_dexpected.md).

## Value

A named list of numeric vectors, one per distribution parameter.
