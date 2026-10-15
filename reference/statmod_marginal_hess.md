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
  ctx = NULL,
  inner = NULL,
  gam = NULL
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

- ctx:

  The evaluation context, or `NULL`.

- inner:

  The inner optimizer, read by the differenced route.

- gam:

  The coefficients the criterion estimates, from
  [`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
  or `NULL`. Their rows and columns follow those of the hyperparameters.

## Value

A square matrix with one row for each row of `idx`, followed by one row
for each coefficient in `gam`, or `NULL` where the penalized information
or the matrix of the determinant cannot be formed or, on a differenced
route, where the two steps disagree.

## Details

The pieces are those of the gradient differentiated once more: the
hyperparameter Hessian of the penalty, the movement of the mode read
through the penalized curvature, and two contributions from the
determinant, one from \\M\\ moving and one from \\K_m\\ moving.

**Contributions of the penalty.** Every contribution of the penalty
comes from a generic of penalties7:
[`penalties7::penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html),
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)
and
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html),
with
[`penalties7::penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.html)
and its second-order companions for a penalty whose Hessian depends on
the coefficients. The assembly therefore contains no branch for a
particular penalty class (a ridge, a random effect, a structured prior).

**Onto the free scale** the chain rule is second order and diagonal,
each hyperparameter having its own link: with \\\theta = h(\eta)\\,
\$\$\partial^2 V/\partial\eta_m\partial\eta_l = h_m'h_l'\\\partial^2
V/\partial\theta_m\partial\theta_l + \delta\_{ml}\\h_m''\\\partial
V/\partial\theta_m.\$\$

**Routes.** For a model with a structural term the function returns
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
where the term implements
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
and a penalty covers the term's own parameters, and
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
otherwise, including the case where that assembly returns no finite
matrix. On the expected information the assembly runs only where
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
admits order 2, and the stencil is used elsewhere. Any other model is
answered by the assembly described above.

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
