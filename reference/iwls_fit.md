# Fit the Smooth Block by Iterated Weighted Least Squares

Runs the scoring iteration on the objective of
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
at fixed hyperparameters, returning the coefficients and the history of
the run.

## Usage

``` r
iwls_fit(obj, start, method, n, pieces_at, verbose = FALSE, groups = NULL)
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

  `TRUE` to print one line per iteration: the objective, the score, the
  step length and the route taken.

## Value

A list of six:

- `par`:

  the stacked coefficients reached, a numeric vector.

- `value`:

  the penalized objective there, unaveraged.

- `converged`:

  a single logical: whether a stopping rule fired, the dimensionless
  verdict included.

- `iterations`:

  how many steps were taken. `1` means the rule was met at the first
  check with no step taken, as a warm start already at the mode reports.

- `score`:

  the score per observation at the point reached.

- `history`:

  a data frame with one row per iteration.

## Two departures from a plain Newton iteration

A curvature that is not positive definite is repaired by flooring its
eigenvalues, and the run continues.
[`solve()`](https://rdrr.io/r/base/solve.html) would signal an error
there and abandon a start that is merely far from the optimum.

The stopping rule is read at the iterate, on a score divided by the
sample size, so a threshold means the same at \\n = 10\\ and at \\n =
10^7\\. The objective itself stays unaveraged, that being the scale the
penalty is added on.

## The dimensionless final verdict

The absolute score of a location equation carries the units \\1/y\\. On
a response three decades small its rounding floor sits above any fixed
threshold, and a run stalled exactly at the optimum reads as a failure.

The verdict therefore adds a second, dimensionless reading, \$\$\max_j
\lvert g_j \rvert / (n\\s_p), \qquad s_p = \sqrt{\mathrm{median}\_j
H\_{jj} / n},\$\$ with \\s_p\\ taken over the equation each coordinate
belongs to.

It can relabel a run that has already stopped as converged, and it can
never stop a run that is still moving. Driving the loop with it was
tried and made the tolerance unreachable at the other end of the scale,
where the objective's magnitude grows with \\\log y\\ and the stall
guard fires first.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the settings,
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
for one step,
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
for the dimensionless reading.
