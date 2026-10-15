# The Exact Outer Hessian of a Model Carrying a Structural Term

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
over the joint vector of coefficients and a filter's own parameters,
which is what the determinant spans there.

## Usage

``` r
statmod_structural_hess(spec, design, coef, hyper, method, idx, basis = NULL)
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

A square matrix on the free scale with one row for each row of `idx`, or
`NULL` where the joint matrix cannot be formed or factorized or the
fifth derivative of the family is unavailable.

## Orders of differentiation

Each order of differentiating the predictor through the recursion
requires one more order of the response's family, because the score that
drives the recursion is read at the predictor that the recursion
produces. The gradient reads \\\partial^3 e/\partial u^3\\ in one
direction through
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html).
The second derivative of the criterion reads \\\partial^4 e/\partial
u^4\\ in two directions through
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
and the fifth derivative of the family, from
[`distributions7::distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.html),
enters with it.

## Form of the assembly

The assembly has the form of
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md).
With \\u\\ the joint vector, \\K\\ the penalized information over it and
\\M\\ the matrix against which the trace is taken,

\$\$\frac{\partial^2 V}{\partial t_m\partial t_l} = -\rho\_{ml} + \hat
b_m^\top K \hat b_l + \tfrac{1}{2}\mathrm{tr}(MK_lMK_m) -
\tfrac{1}{2}\mathrm{tr}\Big(M\frac{\partial K_m}{\partial t_l}\Big),\$\$

with \\\hat b_m = -K^{-1}c_m\\ the mode's movement, \\K_m = S_m +
\partial K/\partial u\[\hat b_m\]\\ and the last trace carrying the
penalty's second derivative, the twice-contracted second derivative of
\\K\\ and the once-contracted one at the second movement of the mode.
Every piece is the one that the coefficient-space assembly uses, read on
the joint vector.

## Blocks that move with their coefficients

A block that moves with its coefficients (`nl()`, `seg()`) beside the
filter contributes nothing here, as it contributes nothing to
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md):
the correction for such a block is written in the coefficient-space
assembly and has no counterpart on the joint vector. Such a model is
admitted at both orders, and the result has the approximation of the
gradient.

## Cost

The assembly makes one call of
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
together with calls of the lower-order recursions, for each pair of
hyperparameters. The stencil refits the model four times for each
hyperparameter.

## See also

[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
for the route that it replaces.
