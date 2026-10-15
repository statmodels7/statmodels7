# The Curvature the Mode Actually Moves By

What separates the true Hessian of the penalized log-likelihood from the
Gauss-Newton matrix that
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
returns, where a term's block moves with its coefficients.

## Usage

``` r
mode_curvature(spec, design, coef, params, npar, offs, total)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at `coef`.

- coef:

  The coefficients at the penalized mode.

- params, npar, offs, total:

  The block bookkeeping.

## Value

A symmetric matrix over the stacked coefficients, with nonzero entries
only on the blocks of refreshable terms.

## Details

\\v = \partial\hat\beta/\partial\theta\\ solves
\\(\partial^2\rho/\partial\beta^2 - \partial^2\ell/\partial\beta^2)v =
-\partial^2\rho/\partial\beta\partial\theta\\, and
\\\partial^2\ell/\partial\beta^2\\ is the true second derivative:
\$\$\sum_i\sum\_{a,b}\ell\_{ab}\frac{\partial\eta_a}{\partial\beta}
\frac{\partial\eta_b}{\partial\beta} + \sum_i\sum_a
\ell_a\frac{\partial^2\eta_a}{\partial\beta^2}.\$\$ The first sum is
what the design gives; the second is zero for every fixed block and is
not for a refreshable one, where \\X\\ is the Jacobian and so
\\\partial^2\eta_i/\partial\beta_c\partial\beta_d = \partial
X\_{id}/\partial\beta_c\\. It is supported on the term's own block, and
that is what keeps this cheap: one call of
[`modelterms7::term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html)
per column, weighted by the score where
[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
weights by \\M\\.

This is the matrix by which the mode moves, and it differs from the
criterion's matrix. The determinant is of the matrix that
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
assembles, and the derivative of the determinant is read from that
matrix; how the mode moves depends on the penalized likelihood and is
read from this one. The two coincide for the observed information when
no block moves with its coefficients. With the expected information a
single matrix for both would give a systematic error in the gradient,
and the two are kept apart.

## See also

[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
[`modelterms7::term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html)
