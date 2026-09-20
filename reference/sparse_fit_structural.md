# Fit a Kinked Penalty Over a Structural Term's Own Parameters

A proximal gradient iteration over the parameters of a score-driven
filter that a penalty with a kink covers, the coefficients and the
term's other parameters held where the rest of the pass left them.

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

Neither route
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
takes for a block of coefficients can be taken here. A coordinate
descent reads the block's columns and its running residual, and a
structural term contributes no column; the proximal branch reads
`beta[block$index]`, and such a penalty covers no entry of the stacked
vector. The parameters are in the design's structural state instead, so
the iteration is written on that subvector: the smooth part is the
negative log-likelihood with this penalty taken back out, its gradient
is
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)'s
with the same correction, and the operator is the penalty's own.

Measured on a panel of eight groups of forty where three carry a level,
the route selects: the coordinates at exactly zero are 0, 0, 1, 3 and 8
of 8 at `lambda` of 0.1, 1, 5, 20 and 80, and the survivors are the
three that carry one. What says the point is the model's is the
Karush-Kuhn-Tucker condition, computed by differencing the smooth
objective with numDeriv and sharing no arithmetic with either the
filter's adjoint or the operator: the stationarity \\\|g +
\lambda\\\mathrm{sign}(b)\|\\ on the coordinates away from zero is
2.1e-08, 2.3e-08 and 6.9e-07 at `lambda` of 1, 5 and 20, and the ones at
zero sit inside the interval the kink opens, 0.83 against 5 and 17.6
against 20.

The hyperparameter is HELD. A criterion that would select it is refused
by
[`assert_criterion_reach()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_reach.md),
a prediction-error criterion reading the degrees of freedom over the
coefficients alone, where such a penalty covers nothing.

THE FIT REPORTS `converged = FALSE` AT A POINT THAT IS AT ITS MODE, and
the flag is not this route's. Traced on the panel above at `lambda` of
20, the passes are: the joint block converges in 58 iterations, this
block runs its budget of 500 and moves the objective by 8.6e-03, then at
the second pass the joint block takes ONE iteration and moves it by
exactly 0 while this block converges in 189. What reports `FALSE` is the
joint block at a pass where there is nothing left to move, which is the
stall guard firing at the mode – the reading
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md)
replaces for availability and which the alternation's own flag does not
yet use. The control is that the same model under
`ridge(~ id, lambda = 20)` and with the levels unpenalized reports
`TRUE`, so the flag follows the extra pass a kinked block costs rather
than anything about the point.

## See also

[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md),
[`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md),
which fits the term's other parameters in the same pass.
