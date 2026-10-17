# Fit a Kinked Penalty Over a Structural Term's Own Parameters

Runs a proximal gradient iteration over the parameters of a score-driven
filter that a penalty with a kink covers, while the coefficients and the
term's other parameters stay where the rest of the pass left them.

## Usage

``` r
sparse_fit_structural(
  obj,
  beta,
  block,
  hyper,
  spec,
  design,
  maxit = 500,
  tol = 1e-08,
  verbose = FALSE
)
```

## Arguments

- obj:

  The full objective.

- beta:

  The stacked coefficients, returned unchanged: this route moves the
  structural state and no coefficient.

- block:

  One entry of `statmod_blocks()$sparse`, whose `structural` is `TRUE`.

- hyper:

  The hyperparameters.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose structural state this writes into.

- maxit:

  The iteration budget.

- tol:

  The stopping tolerance.

- verbose:

  Whether the optimizer prints its own trace.

## Value

A list shaped like
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)'s.

## Details

The two routes that
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
takes for a block of coefficients cannot be taken here. A coordinate
descent reads the block's columns and its running residual, and a
structural term contributes no column; the proximal branch reads
`beta[block$index]`, and such a penalty covers no entry of the stacked
vector. The parameters are in the design's structural state instead, so
the iteration is written on that subvector: the smooth part is the
negative log-likelihood from which this penalty is taken back out, its
gradient is
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)'s
with the same correction, and the operator is the penalty's own.

The iteration stops when the gradient criterion `crit_grad(tol)` is met.
A coordinate at exactly zero is one whose gradient of the smooth part
lies inside the interval that the kink opens.

The hyperparameter is held. A prediction-error or cross-validation
criterion
([`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md))
that would select it is rejected by
[`assert_criterion_reach()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_reach.md),
because such a criterion works only over the coefficients, where this
penalty covers nothing.

The convergence flag of a fit that uses this route is the flag of the
alternation. Beside a structural term or a kinked block the alternation
requires the point to be settled in the sense of
[`alternation_settled()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md)
(the Newton decrement of the smooth part and the decrement of the kinked
zeros both at most
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)).

## See also

[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md),
[`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md),
which fits the term's other parameters in the same pass.
