# Estimate the Hyperparameters

Runs `outer_optimizer` on the marginal criterion, refitting the
coefficients at every hyperparameter it tries.

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
  vb
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

  An optimizers7 optimizer, or `NULL` to let the availability of the
  exact gradient decide.

- beta:

  The starting coefficients, stacked.

- approx:

  The approximation for the expected information.

- maxit, tol:

  The alternation's budget and tolerance.

- vb:

  The resolved verbosity.

## Value

A list with `par`, `hyper`, `value`, `criterion`, `converged`, `history`
and the inner results.

## Details

Each evaluation is a whole inner fit, warm-started from the previous
one, and that is what makes the search affordable: after the first few
hyperparameters the coefficients move very little and the inner loop
converges in two or three iterations.

**The optimizer is chosen by whether the gradient exists.** Where
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
applies – the observed information, and penalties whose Hessian is
linear in their hyperparameters – the criterion is handed its exact
derivative and
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
is the default; otherwise the search compares values and
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html)
is. An optimizer given explicitly is used as given, and one that needs a
gradient it cannot be given will say so itself.

## See also

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
