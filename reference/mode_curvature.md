# The Curvature the Mode Actually Moves By

What separates the true Hessian of the penalized log-likelihood from the
Gauss-Newton matrix
[`statmod_information_at`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
returns, where a term's block moves with its coefficients.

## Usage

``` r
mode_curvature(spec, design, coef, params, npar, offs, total)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at `coef`.

- coef:

  The coefficients at the penalized mode.

- params, npar, offs, total:

  The block bookkeeping.

## Value

A square matrix, zero everywhere no refreshable term reaches.

## Details

\\v = \partial\hat\beta/\partial\theta\\ solves
\\(\partial^2\rho/\partial\beta^2 - \partial^2\ell/\partial\beta^2)v =
\partial^2\rho/\partial\beta\partial\theta\\, and
\\\partial^2\ell/\partial\beta^2\\ is the TRUE second derivative:
\$\$\sum_i\sum\_{a,b}\ell\_{ab}\frac{\partial\eta_a}{\partial\beta}
\frac{\partial\eta_b}{\partial\beta} + \sum_i\sum_a
\ell_a\frac{\partial^2\eta_a}{\partial\beta^2}.\$\$ The first sum is
what the design gives; the second is zero for every fixed block and is
not for a refreshable one, where \\X\\ is the Jacobian and so
\\\partial^2\eta_i/\partial\beta_c\partial\beta_d = \partial
X\_{id}/\partial\beta_c\\. It is supported on the term's own block,
which is what keeps this cheap: one call of
[`term_block_contract`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html)
per column, weighted by the SCORE where
[`u_refresh`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
weights by \\M\\.

⚠️ It is the mode's matrix and NOT the criterion's. The determinant is
of whatever
[`statmod_information_at`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
assembles, and its derivative reads that one; how the mode MOVES is a
fact about the penalized likelihood and reads this one. Confusing the
two is the defect this file already records for the expected
information, in a second place.

Measured on `nl(a ~ 0 + ridge(~grp))` against a finite difference of the
criterion with the mode refitted: the gradient goes from 1.5e-04 to
8.8e-09 at 320 observations and from 1.6e-05 to 1.6e-09 at 960, on the
observed route and the expected one alike.

## See also

[`u_refresh`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
[`term_block_contract`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html)
