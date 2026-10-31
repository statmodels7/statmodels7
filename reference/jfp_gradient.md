# The Gradient of the Jeffreys Penalty in the Joint Vector

Computes \\\partial_k\sum_b a_b\log\|\tilde I_b\|\\ over every
coordinate of \\u\\.

## Usage

``` r
jfp_gradient(st)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

## Value

A numeric vector over \\u\\.

## Details

For a block whose rows do not move, \$\$\partial_k\log\|\tilde I\| =
\mathrm{tr}(\tilde I^{-1}\partial_k\tilde I) = \sum_t w_t\sum_c
V\_{c,t,k}\\s\_{c,t},\qquad s\_{c,t} =
-\sum\_{q,r}h\_{qr,t}\\\partial\_{\eta_c}E\_{qr,t},\$\$ with
\\V\_{c,t,k}\\ the derivative of the predictor of equation \\c\\ in
\\u_k\\ and \\\partial\_{\eta_c}E\_{qr}\\ from
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html).
This is the contraction of
[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
with the block's leverages in place of the diagonal of a matrix over
\\u\\, which cannot hold the constant direction.

## See also

[`jfp_value()`](https://statmodels7.github.io/statmodels7/reference/jfp_value.md),
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
