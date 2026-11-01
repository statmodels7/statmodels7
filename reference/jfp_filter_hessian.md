# The Hessian's Shares From a Filtered Block's Own Movement

Computes the terms of \\\partial\_{kl}\log\|\tilde I\|\\ that the
movement of the block's rows contributes, for a block inside the
subformula of a filter's parameter, whose rows are derivatives of the
filtered predictors.

## Usage

``` r
jfp_filter_hessian(st, bs, dE)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- bs:

  A block's entry of the state, of kind `"filter"`.

- dE:

  The derivative of the expected information, from
  [`jfp_dexpected()`](https://statmodels7.github.io/statmodels7/reference/jfp_dexpected.md).

## Value

A list with `H`, a matrix over the joint vector, and `BV`, the rows to
add to the low-rank factor, in the convention of
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md).

## Details

With the whitened directions \\d_i = C U^{-1}e_i\\ the block's rows are
\\Z\_{q,t,i} = \partial\eta\_{q,t}\[d_i\]\\, \\C\tilde I^{-1}C^\top =
\sum_i d_id_i^\top\\, and the derivative of the rows in \\u_k\\ is
\\R^{(i)}\_{q,t,k} = \partial^2\eta\_{q,t}\[d_i, k\]\\. Then
\\\mathrm{tr}(\tilde I^{-1}\partial\_{kl}\tilde I)\\ gains \$\$2\sum_t
w_t\sum\_{q,r}\sum_i\Big\[\sum_c\partial_c\mathcal I\_{qr}
\big(V\_{c,t,k}Z\_{q,t,i}R^{(i)}\_{r,t,l} + (k\leftrightarrow l)\big) +
\mathcal I\_{qr}Z\_{q,t,i}\\\partial^3\eta\_{r,t}\[d_i,k,l\] + \mathcal
I\_{qr}R^{(i)}\_{q,t,k}R^{(i)}\_{r,t,l}\Big\],\$\$ and the whitened
derivative \\U^{-\top}\partial_k\tilde I U^{-1}\\ gains
\\X^{(a)}\_{kb} + X^{(b)}\_{ka}\\ with \\X^{(a)}\_{kb} = \sum_t w_t
\sum\_{q,r}\mathcal I\_{qr}R^{(a)}\_{q,t,k}Z\_{r,t,b}\\. One call of
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
per direction gives both \\\partial^3\eta\[d_i,\cdot,\cdot\]\\ and
\\R^{(i)}\\, so the second derivatives of the rows are never held for
more than one direction.
