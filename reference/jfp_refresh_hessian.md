# The Hessian's Shares From a Refreshed Block's Own Movement

Computes the terms of \\\partial\_{kl}\log\|\tilde I\|\\ that the
movement of the block's rows contributes, for a block inside a
subformula of a term whose block is recomputed at the coefficients.

## Usage

``` r
jfp_refresh_hessian(st, bs, dE)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- bs:

  A block's entry of the state, of kind `"refresh"`.

- dE:

  The derivative of the expected information, from
  [`jfp_dexpected()`](https://statmodels7.github.io/statmodels7/reference/jfp_dexpected.md).

## Value

A list with `H`, a matrix over the joint vector, and `BV`, the rows to
add to the low-rank factor at the term's positions.

## Details

With \\\tilde I = \sum_t w_tI_t\tilde V_t^\top\tilde V_t\\ over the
block's equation, \\P = \tilde V\tilde I^{-1}\\ and \\I\_{c,t}\\ the
derivative of \\I_t\\ in the predictor of equation \\c\\,
\\\mathrm{tr}(\tilde I^{-1}\partial\_{kl}\tilde I)\\ gains \$\$2\sum_t
w_t\Big\[\sum_c I\_{c,t}\big(X_c\[t,k\]R\_{tl} +
X_c\[t,l\]R\_{tk}\big) +
I_t\\\mathrm{tr}\big(P_t^\top\\\partial\_{kl}\tilde V_t\big) +
I_t\\\partial_k\tilde V_t\\\tilde I^{-1}\partial_l\tilde
V_t^\top\Big\],\$\$ with \\R\_{tl} = P_t\\\partial_l\tilde V_t^\top\\;
the rows move in the term's own coefficients only.
\\\partial\_{kl}\tilde V\\ is read from
[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html)
by the symmetry of the third derivative of the predictor. The product
\\\mathrm{tr}(\tilde I^{-1}\partial_k\tilde I\\\tilde I^{-1}
\partial_l\tilde I)\\ gains the whitened derivative \\\sum_t
w_tI_t(W\_{k,t}^\top Z_t + Z_t^\top W\_{k,t})\\, \\W_k =
\partial_k\tilde V\\U^{-1}\\, in the rows of the low-rank factor.
