# Estimate the Hyperparameters

Runs `outer_optimizer` on the marginal criterion, refitting the
coefficients at every hyperparameter that it tries.

## Usage

``` r
outer_fit(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  method,
  optimizer,
  beta,
  approx,
  maxit,
  tol,
  vb,
  hold_hyper = FALSE
)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- blocks:

  The block split.

- hyper:

  The starting hyperparameters.

- inner_optimizer:

  How the smooth block is fitted.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- optimizer:

  An optimizers7 optimizer, or `NULL` to choose one with
  [`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md)
  from the derivatives that the criterion supplies. With `NULL`, the
  stopping rule and the line search are also set by the function.

- beta:

  The starting coefficients, stacked.

- approx:

  The approximation for the expected information.

- maxit, tol:

  The alternation's budget and tolerance.

- vb:

  The resolved verbosity.

- hold_hyper:

  `TRUE` to hold the hyperparameters at `hyper` and search over the
  criterion's own coefficients alone.

## Value

A list with `par`, `hyper`, `value`, `criterion`, `converged`,
`settled`, `obj`, `hist_blocks`, `hist_inner`, `aliased`, `spec`,
`hist_outer`, `iterations`, `evaluations`, `exact_gradient`,
`exact_hessian`, `optimizer`, `held_coef`, `marginal` and `pinned`.
`hist_outer` holds the trace of the search.

## Details

Each evaluation is a whole inner fit, warm-started from the incumbent,
the best usable point found so far. Near the optimum the coefficients
move little, so the inner fit converges in few iterations.

The optimizer is chosen by
[`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md)
from the derivatives that the criterion supplies:
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
where the exact Hessian is used,
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
where the gradient is exact (or a covariance class spans a structural
term's parameters),
[`optimizers7::brent()`](https://statmodels7.github.io/optimizers7/reference/brent.html)
for a single coordinate without a gradient, and
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html)
otherwise. An optimizer given explicitly is used as given, and an
optimizer that needs a gradient that is not available signals its own
error.

## See also

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
