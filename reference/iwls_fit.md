# Fit the Smooth Block by Iterated Weighted Least Squares

Runs the scoring iteration on the objective of
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
at fixed hyperparameters, returning the coefficients and the history of
the run.

## Usage

``` r
iwls_fit(
  obj,
  start,
  method,
  n,
  pieces_at,
  verbose = FALSE,
  groups = NULL,
  frozen = integer(0),
  backup_at = NULL,
  damp_on_reject = TRUE,
  kinks_at = NULL,
  newton_at = NULL,
  switch_at = NULL,
  switch_after = Inf
)
```

## Arguments

- obj:

  The objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it.

- start:

  The starting coefficients, stacked into one numeric vector.

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object, carrying the curvature, the decomposition, the budget and the
  stopping rule.

- n:

  The number of observations, the divisor of the scaled stopping rule.

- pieces_at:

  A function of the stacked coefficients returning what
  [`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  needs, as
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds it.

- verbose:

  `TRUE` to print one line per iteration: the objective, the score per
  observation and the step length.

- groups:

  A list of integer index vectors, one per equation, in the numbering of
  the coefficients handed in. `NULL` treats them as one group, which is
  right where the whole vector is being fitted.

- frozen:

  Positions held at their starting values, an integer vector in the same
  numbering. They are dropped from every solve, so the run maximizes
  over the rest, and the stopping rule reads only the free coordinates,
  because the constrained optimum leaves in a held coordinate a score
  that does not vanish.

- backup_at:

  `NULL`, or a function, like `pieces_at`, that builds the pieces on the
  expected information. The expected information takes the step where
  the observed step finds no acceptable point. An observed penalized
  information that is not positive definite is first repaired by
  [`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
  and its step tried. A trial point whose objective raises an error is
  read as rejected.
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  passes one for a method that
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  settled with the fallback.

- damp_on_reject:

  Whether a line search that finds no acceptable step raises the
  Levenberg damping and retries (the default) or ends the run.
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  passes `FALSE` where the design carries a sharp break-point term,
  whose objective has a kink in the break-point at the observations:
  there every Gauss-Newton step from the minimum is rejected, and
  escalating the damping at every inner fit costs iterations without
  reaching anything that a stop would not.

- kinks_at:

  `NULL`, or a function of the coefficients returning the rows of the
  linear constraints that hold the objective on its kinks there, a
  matrix with one column per coefficient, as
  [`kink_constraints()`](https://statmodels7.github.io/statmodels7/reference/kink_constraints.md)
  finds them. Where a line search rejects every step length, the problem
  is solved on the subspace that those constraints and `frozen` leave
  free
  ([`iwls_tied()`](https://statmodels7.github.io/statmodels7/reference/iwls_tied.md));
  from the point reached one step on every coordinate is tried, and
  where it is rejected too the run stops there, with the verdict of the
  constrained fit.

- newton_at:

  `NULL`, or a function of the coefficients returning the full Hessian
  of the objective, the scoring pieces' curvature plus the term that a
  block moving with its coefficients adds. Where the scoring step is
  rejected or shrunk below a tenth, the Newton step on it is tried by
  [`iwls_newton_step()`](https://statmodels7.github.io/statmodels7/reference/iwls_newton_step.md)
  and the one that decreases the objective more is taken.

- switch_at:

  `NULL`, or a function, like `pieces_at`, that builds the pieces on the
  observed information. After `switch_after` iterations without
  convergence, or earlier when
  [`iwls_slow()`](https://statmodels7.github.io/statmodels7/reference/iwls_slow.md)
  reports a slow contraction of the score, it takes the place of
  `pieces_at`, and `pieces_at` that of `backup_at`.
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  passes one for a method that
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  settled on the expected information from `"auto"`.

- switch_after:

  The iteration after which `switch_at` replaces `pieces_at`.

## Value

A list of ten: the six below; `note`, the reason for which a run stopped
or `NULL`; `aliased`, the coordinates that the pivot left out;
`switched`, the iteration at which the run moved to `switch_at`, or 0;
and `fallback`, a named integer vector of four counts. The count
`indefinite` is the number of iterations at which the observed penalized
information was not positive definite, `search` the number at which the
expected pieces stepped in for a step that found no acceptable point,
`error` the number of trial points whose objective raised, and `newton`
the number of iterations at which the Newton step on the full Hessian
was taken in place of the scoring step. The six are:

- `par`:

  the stacked coefficients reached, a numeric vector.

- `value`:

  the penalized objective there, not averaged over the observations.

- `converged`:

  a single logical: whether a stopping rule fired, including the
  dimensionless verdict.

- `iterations`:

  the index of the iteration at which the loop ended. `1` means the rule
  was met at the first check, as a warm start already at the mode
  reports; the full step taken after the rule is met is not counted.

- `score`:

  the score per observation at the point reached.

- `history`:

  a data frame with one row per accepted step, or `NULL` when no step
  was accepted.

## Repair of the curvature

A curvature that is not positive definite is repaired by flooring its
eigenvalues, and the run continues.
[`solve()`](https://rdrr.io/r/base/solve.html) would signal an error
there and abandon a start that is merely far from the optimum.

## The stopping rule

The built-in stopping rule compares a score divided by the sample size
with `tol`, so that a threshold has the same meaning at \\n = 10\\ and
at \\n = 10^7\\, and it also requires the Newton decrement to be at most
\\10^{-6}\\. After the rule is met, one more full step is taken. The
objective itself is not averaged, because the penalty is added on that
scale.

## The dimensionless final verdict

The absolute score of a location equation carries the units \\1/y\\. On
a response three decades small its rounding floor sits above any fixed
threshold, and a run stalled exactly at the optimum reads as a failure.

The verdict therefore includes a second, dimensionless reading,
\$\$\max_j \lvert g_j \rvert / (n\\s_p), \qquad s_p =
\sqrt{\mathrm{median}\_j H\_{jj} / n},\$\$ where \\s_p\\ is taken over
the equation to which each coordinate belongs.

It can relabel a run that has already stopped as converged, and it never
stops a run that is still moving. Driving the loop with it would make
the tolerance unreachable at the other end of the scale, where the
magnitude of the objective grows with \\\log y\\ and the stall guard
fires first.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the settings,
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
for one step,
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
for the dimensionless reading.
