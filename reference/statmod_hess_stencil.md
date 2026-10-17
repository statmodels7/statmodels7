# The Outer Hessian by One Difference of the Exact Gradient

Returns the criterion's second derivative in the hyperparameters where
there is no analytic route: one central difference of
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
in which the coefficients are refitted at every probe.

## Usage

``` r
statmod_hess_stencil(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  inner = NULL,
  h = hess_stencil_step(),
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

  The coefficients from which the probes start.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The hyperparameter index, from
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- basis:

  The integrated basis, or `NULL`.

- inner:

  The inner optimizer with which the probes refit;
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  where the argument is `NULL`.

- h:

  The step, on the free scale on which the search runs. For a
  coefficient that the criterion estimates the step is \\h \max(1,
  \|\gamma_j\|)\\.

- gam:

  The coefficients that the criterion estimates, from
  [`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
  or `NULL`. Each probe holds them at its own values and refits the
  rest, as
  [`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
  does.

## Value

The Hessian on the free scale over the hyperparameters followed by the
coefficients in `gam`, or `NULL` where a probe could not be evaluated or
the two steps disagree.

## Details

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
uses this function for a model with a structural term for which the
joint assembly
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
does not apply: a filter without a penalty over its own parameters, a
term of the likelihood shape such as `regime()`, or a term that does not
implement
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html).
It also uses it on the expected information where
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
does not admit order 2. The analytic assembly over the coefficients
alone is not the criterion's Hessian when a filter's own parameters move
with the hyperparameter.
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
covers that case with the fourth derivative of the recursion from
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
and the fifth derivative of the family from
[`distributions7::distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.html).

The function differences the exact gradient once, so the result contains
a single layer of finite differences. The matrix is symmetrized by
averaging it with its transpose.

## Probes

Each probe refits the coefficients from the coefficients supplied and
from the structural state as it stands, which is restored before each
probe, so that the two probes of a column differ only in the
hyperparameter. The location error of the mode is nearly the same at
\\+h\\ and \\-h\\ and cancels in the difference. Where the inner
optimizer is
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
the probes use its tolerance `min(tol, hess_stencil_inner_tol())`, a
budget of at least 500 iterations and no stopping criterion of its own;
any other inner optimizer is used as given.

## Rejection of an unresolved stencil

A stencil inherits the reproducibility of the quantity that it
differences, divided by the step, so where the mode is poorly located
the result is noise. The matrix is therefore computed at `h` and at
`3 * h`, and the function returns `NULL` where the largest absolute
difference between the two exceeds
[`hess_stencil_tol()`](https://statmodels7.github.io/statmodels7/reference/hess_stencil_step.md)
times the largest absolute entry of the first. Callers treat a `NULL`
result as an unavailable curvature.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`hess_stencil_step()`](https://statmodels7.github.io/statmodels7/reference/hess_stencil_step.md)
