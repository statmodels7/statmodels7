# The Outer Criterion on the Compiled Route

`engine_outer_at()` returns the compiled context of an outer point, from
which
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
and
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
read the criterion, its gradient and its Hessian on the models that the
compiled route covers, or `NULL` where the R route computes them.

## Usage

``` r
engine_outer_at(ctx, spec, design, coef, hyper, method, basis = NULL)

engine_outer_build(spec, design, coef, hyper, basis)

engine_outer_pen(spec, design, coef, hyper)

engine_outer_members(spec, design, idx, hyper)

engine_ctx_solve(ctx, x)

engine_flags

engine_outer_complete(eo, spec, design, coef, hyper)

engine_outer_mvt(u, un, coef, spec, hyper)

engine_outer_chart(u, hyper)

outer_search_settings(optimizer, nv)

engine_pen_spec(spec, design, idx, hyper, beta, how)

engine_chart_static(s, inverted = FALSE)

engine_mvt_static(u, ul)
```

## Format

`engine_flags`: an environment of switches of the compiled route.

## Arguments

- ctx:

  An evaluation context from
  [`outer_context()`](https://statmodels7.github.io/statmodels7/reference/outer_context.md),
  or `NULL`.

- spec, design, coef, hyper, method, basis:

  As for
  [`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md).

- idx:

  The outer index.

- x:

  A vector of the length of the stacked coefficients.

- eo:

  An outer context of `engine_outer_at()`.

- u:

  A penalized unit of
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  whose penalty is a multivariate gaussian prior.

- un:

  The unit's description by
  [`engine_penalty_units()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md).

- optimizer:

  An optimizer of optimizers7.

- nv:

  The length of the outer vector.

- beta:

  The stacked coefficients at the start.

- how:

  The solver of the inner fit, as for
  [`engine_model()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md).

- s:

  A matrix parameter of parameters7 (a covariance chart).

- inverted:

  Whether the chart is of the precision.

- ul:

  The unit's description by
  [`engine_penalty_units()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md).

## Value

`engine_outer_at()`: an external pointer, or `NULL`.
`engine_outer_members()`: a list with `mrow`, `munit`, `mchain`, `h1`
and `h2`, or `NULL` where a member is not covered.

## Details

A point is covered when `spec@engine` is `"auto"`, the method is
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
on the observed information, the model is one that
[`engine_model()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md)
describes, the family carries its derivatives to the fourth order in the
registry, and every penalized unit is a quadratic penalty or a gaussian
prior with a fixed mean and no map. The compiled context computes, at
the mode, the derivatives of the log-density to the fourth order on the
link scale, the observed information, the penalized information pinned
at
[`pinned_coords()`](https://statmodels7.github.io/statmodels7/reference/pinned_coords.md)
with its inverse, the matrix of the traces and the leverages, as the
evaluation context of
[`outer_context()`](https://statmodels7.github.io/statmodels7/reference/outer_context.md)
holds them. Where the R route would take a branch that the compiled
context does not reproduce (a curvature that is not finite, a Cholesky
factor that fails or whose condition is within a factor of 100 of the
acceptance threshold of
[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md),
and therefore the eigendecomposition and the retry over
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)),
the context reports it and the R route computes the point. The context
is stored in `ctx`, keyed by the integrated basis, so that the
criterion, the gradient and the Hessian at one point read one context.

`engine_outer_members()` describes the outer index to the compiled
context: for each member line of
[`index_members()`](https://statmodels7.github.io/statmodels7/reference/index_members.md)
its row, the position of its penalized unit in
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
and the derivative of its hyperparameter's inverse link at its value;
and, for each row, the first and second derivatives of the inverse link
at the first member's value, which the Hessian reads.

`engine_outer_pen()` returns the model and the penalty as the outer
context reads them, `list(model, pen)`, or `NULL` where the point is not
covered:
[`engine_model()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md)
on the Cholesky route, with the pieces that the context reads beside the
penalty's units (the positions of a t prior's hyperparameters, an
additive penalty's blocks, a covariance chart).

`engine_ctx_solve()` returns \\K^{-1}x\\ with \\K\\ the mode's pinned
penalized information, from the compiled context that `ctx` already
holds, or `NULL` where `ctx` holds none or the context does not cover
the point, in which case the caller reads the inverse of
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md).
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md)
and
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
read it. The switch `engine_flags$mode_inverse` turns it off.

`engine_outer_complete()` gives the context the derivatives of each
multivariate t prior in its hyperparameters (`engine_outer_mvt()`),
which the criterion does not read, before the first gradient or Hessian.
It returns the context, or `NULL` where a prior is not covered.

`engine_outer_mvt()` describes a multivariate Student t prior to the
compiled context at the mode: the derivatives of the penalty in its
hyperparameters
([`penalties7::penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html),
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)),
the blocks of
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html)
and
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html)
group by group, whether
[`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html)
holds, and per hyperparameter the derivative of the inverse scale matrix
(zero for `nu`) and of `nu + p` (one for `nu`), from which the compiled
context forms the movement of the penalty's Hessian with the
coefficients. It returns `NULL` where the chart's free values do not
match the hyperparameters.

`engine_outer_chart()` describes the covariance chart of a multivariate
gaussian prior to the compiled context: a code (1 for
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html),
of the covariance or of the precision, 2 for
[`parameters7::compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.html),
3 for
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)),
the free values of the chart in its own order, which is the order of the
penalty's hyperparameters, and, for the last two, the inverse links of
the scale and of the correlation with their first two derivatives at the
free values. It returns `NULL` for any other chart.

`outer_search_settings()` returns the settings of a descent method for
the compiled outer search
([`optimizers7::descent_settings()`](https://statmodels7.github.io/optimizers7/reference/descent_settings.html)),
or `NULL` where the search runs in R: an optimizer that does not run the
descent loop, and a stopping rule without a compiled form. The switch
`engine_flags$outer_search` turns the compiled search on,
`engine_flags$outer_point` the point computed in compiled code within
it, and `engine_flags$outer_pen` the penalty built in compiled code at
each point; all three are on.

`engine_pen_spec()` describes the penalty to the compiled construction
of the outer search (`engine_pen_build()`): per penalized unit its kind,
its positions, its description by
[`engine_penalty_units()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md)
for the inner fit and for the outer context, its hyperparameters in the
order of the penalty's parameters, and the static pieces of its kind;
and the members of the outer index with the links of its rows. It
returns `NULL` where a unit is not one that the compiled construction
writes (kinds 1 to 4, a kind 4 unit on a log-Cholesky, compound symmetry
or AR(1) chart, no map on kinds 3 and 4), or where a link has no
compiled route.

`engine_chart_static()` describes a covariance chart to the compiled
construction (`engine_charts.h`): its code (1 log-Cholesky, 2 compound
symmetry, 3 AR(1), 4 `dr_prod()` over a `correlation_matrix()`), its
dimension and number of free values, the orientation, the positions of
the log-Cholesky factor or of the angles (0-based), and the links of the
scales and of the correlation. It returns `NULL` for any other chart,
and for a `dr_prod()` whose links are not the log and the doubly bounded
link. `engine_mvt_static()` describes a multivariate t prior: its chart,
its location, the role of each hyperparameter (the chart's free value,
0-based, or -1 for `nu`), and for an odd dimension the coefficients of
the series of `mvt_cp_d()`.

## See also

[`engine_inner_fit()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md),
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
