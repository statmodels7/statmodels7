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
  inner_method,
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

- inner_method:

  How the smooth block is fitted.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- optimizer:

  An optimizers7 optimizer.

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
one, which is what makes the search affordable: after the first few
hyperparameters the coefficients move very little and the inner loop
converges in two or three iterations.

The default optimizer compares values and asks for no derivative. An
exact gradient of this criterion is available in principle – the
envelope theorem kills the term through \\d\hat\beta/d\theta\\ in the
first two pieces, so only the determinant's derivative is left – and it
is not written; a gradient-based optimizer supplied here is handed a
numerical one, whose step has to stay well above the noise the inner
tolerance leaves on the criterion.

## See also

[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
