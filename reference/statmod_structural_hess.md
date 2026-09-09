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

A square matrix on the free scale, one row per row of `idx`, or `NULL`
where the joint matrix could not be formed.

## Why it needs a fourth order

Each order of differentiating the predictor through the recursion pulls
in one more order of the response's family, the score the recursion is
driven by being read at the predictor it produces. The gradient reads
\\\partial^3 e/\partial u^3\\ in one direction through
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html);
the criterion's own second derivative reads \\\partial^4 e/\partial
u^4\\ in two, through
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
and the family's FIFTH derivative with it.

## The shape is [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)'s

With \\u\\ the joint vector, \\K\\ the penalized information over it and
\\M\\ the matrix the trace is taken against,

\$\$\frac{\partial^2 V}{\partial t_m\partial t_l} = -\rho\_{ml} + \hat
b_m^\top K \hat b_l + \tfrac{1}{2}\mathrm{tr}(MK_lMK_m) -
\tfrac{1}{2}\mathrm{tr}\Big(M\frac{\partial K_m}{\partial t_l}\Big),\$\$

with \\\hat b_m = -K^{-1}c_m\\ the mode's movement, \\K_m = S_m +
\partial K/\partial u\[\hat b_m\]\\ and the last trace carrying the
penalty's second derivative, the twice-contracted second derivative of
\\K\\ and the once-contracted one at the mode's second movement. Every
piece is the one the coefficient-space assembly uses, read on the joint
vector.

## What it does not carry

A block that MOVES with its coefficients – `nl()`, `seg()` – beside the
filter contributes nothing here, exactly as it contributes nothing to
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md):
that correction is written in the coefficient-space assembly and has no
joint twin. Such a model is admitted at both orders and the
approximation is the gradient's own.

## Cost

One
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
and three lower recursions per PAIR of hyperparameters, against the
stencil's four refits per hyperparameter.

## See also

[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
for the route it replaces.
