# The Hessian of the Jeffreys Penalty in the Joint Vector

Computes \\\partial\_{kl}\sum_b a_b\log\|\tilde I_b\|\\ as a matrix with
the sparsity of the likelihood's Hessian minus a product of low rank.

## Usage

``` r
jfp_hessian(st, assemble = TRUE)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- assemble:

  Whether to return the sum as a dense matrix as well.

## Value

A list with `sparse` (the first term times the powers, in the storage
type of the information), `lowrank` (\\B\\ times the square roots of the
powers, one block of columns per Jeffreys block), `weights` (the
per-observation weights of the first term, keyed as
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
reads its `H`, so that `sparse` is that function's result at
`H = weights` plus `refresh`), `refresh` (the part that a refreshed
block contributes, or `NULL`) and, when `assemble` is `TRUE`, `hessian`,
equal to `sparse - tcrossprod(lowrank)` as a base matrix.

## Details

For a block whose rows do not move, \$\$\partial\_{kl}\log\|\tilde I\| =
\mathrm{tr}(\tilde I^{-1} \partial\_{kl}\tilde I) - \mathrm{tr}(\tilde
I^{-1}\partial_k\tilde I\\ \tilde I^{-1}\partial_l\tilde I).\$\$ The
first term is \\\sum\_{c,d}X_c^\top\mathrm{diag}(w\\W\_{cd})X_d\\ with
\\W\_{cd} = -\sum\_{q,r}h\_{qr}\\\partial^2\_{\eta_c\eta_d}E\_{qr}\\,
from
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html),
assembled by
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md).
Where an equation carries a term whose block is recomputed at the
coefficients, the second derivative of its predictor adds \\\sum_t
w_t\sum_c s\_{c,t}\\\partial^2\eta\_{c,t}/\partial u_k \partial u_l\\,
with \\s_c\\ the weights of
[`jfp_gradient()`](https://statmodels7.github.io/statmodels7/reference/jfp_gradient.md),
read through
[`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md).

The second term is \\BB^\top\\: row \\k\\ of \\B\\ holds the upper
triangle of \\U^{-\top}\partial_k\tilde I\\U^{-1}\\, the off-diagonal
entries multiplied by \\\sqrt2\\, so that its rank is at most \\\tilde
p(\tilde p+1)/2\\.

## See also

[`jfp_value()`](https://statmodels7.github.io/statmodels7/reference/jfp_value.md),
[`jfp_gradient()`](https://statmodels7.github.io/statmodels7/reference/jfp_gradient.md)
