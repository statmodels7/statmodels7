# The Exact Hessian of the Marginal Criterion

\\\partial^2 V/\partial\eta^2\\ at the penalized mode, over the free
scale of the hyperparameters under estimation.

## Usage

``` r
statmod_marginal_hess(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  ctx = NULL
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

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The outer index.

- basis:

  The integrated subspace, or `NULL`.

## Value

A square matrix, one row per row of `idx`, or `NULL` where the
determinant does not exist.

## Details

The pieces are those of the gradient differentiated once more: the
hyperparameter Hessian of the penalty, the movement of the mode read
through the penalized curvature, and two contributions from the
determinant, one from \\M\\ moving and one from \\K_m\\ moving.

**Everything the penalty contributes is asked of the penalty**
([`penalties7::penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html),
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html),
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)),
so a penalty that is not quadratic in its hyperparameters is covered by
the same assembly with no branch here: a ridge, a random effect, a
structured prior.

**Onto the free scale** the chain rule is second order and diagonal,
each hyperparameter having its own link: with \\\theta = h(\eta)\\,
\$\$\partial^2 V/\partial\eta_m\partial\eta_l = h_m'h_l'\\\partial^2
V/\partial\theta_m\partial\theta_l + \delta\_{ml}\\h_m''\\\partial
V/\partial\theta_m.\$\$

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
