# How a Moving Penalty Hessian Enters the Criterion's Second Derivative

Computes the part of \\\partial K_m/\partial t_l\\ that a penalty whose
Hessian depends on the coefficients contributes, placed in the stacked
coefficient space, where \\K_m = S_m + T\[b_m\]\\ is the determinant's
matrix differentiated in hyperparameter \\m\\.

## Usage

``` r
statmod_penalty_second(
  spec,
  design,
  coef,
  hyper,
  idx,
  m,
  l,
  bm,
  bl,
  bml,
  total,
  M = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- idx:

  The outer index.

- m, l:

  The two rows of the index.

- bm, bl, bml:

  The mode's movement along each and along the pair.

- total:

  The stacked width.

- M:

  `NULL`, or the matrix against which the result is traced; the trace is
  then returned instead of the matrix, and a Jeffreys penalty's share is
  computed without assembling it.

## Value

A `total` square matrix, or `NULL` if no penalty's Hessian moves with
the coefficients; with `M`, a single number.

## Details

Differentiating \\S_m(\beta, t) + \partial S/\partial\beta\\\[b_m\]\\
once more along \\t_l\\, while \\\beta\\ moves by \\b_l\\, gives
\$\$\frac{\partial S_m}{\partial\beta}\[b_l\] + \frac{\partial
S_l}{\partial\beta}\[b_m\] + \frac{\partial^2 S}{\partial\beta^2}\[b_l,
b_m\] + \frac{\partial S}{\partial\beta}\[b\_{ml}\]\$\$ beside
\\S\_{ml}\\, which
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
already carries. The first two are the mixed derivative \\\partial^2
S/\partial\beta\\\partial\theta\\ taken along two directions, \\b_l\\
for the hyperparameters of row \\m\\ and \\b_m\\ for those of row \\l\\;
each is read from
[`penalties7::penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.html)
for the hyperparameters that the penalty owns. The third is
[`penalties7::penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.html)
and the fourth
[`penalties7::penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.html).
A penalty for which
[`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html)
returns `TRUE` does not contribute, and a kinked penalty, a mixed
covariance class and a penalty over a structural term's own parameters
do not contribute either, so for a model carrying only those the
function returns `NULL` and the Hessian is left untouched.

## See also

[`statmod_penalty_dbeta()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_dbeta.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`statmod_pe_derivs()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe_derivs.md)
