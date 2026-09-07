# Refit With One Coefficient Held at a Value

Maximizes the penalized likelihood over every coefficient but one, which
is held at `value`. What comes back is the restricted maximum, the
coefficients that attain it and the score there.

## Usage

``` r
statmod_restrict(fit, param, coefname, value)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The distribution parameter whose equation carries the coefficient, a
  single string.

- coefname:

  The coefficient's name, a single string, as
  [`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
  reports it.

- value:

  The value to hold it at, a single number.

## Value

A list:

- `coefficients`:

  the restricted estimates, a named list.

- `loglik`:

  the log-likelihood there, unpenalized.

- `objective`:

  the penalized objective there.

- `score`:

  the penalized score at the restricted point, a numeric vector over the
  stacked coefficients. Its free entries are zero to the tolerance; the
  held one is what the statistics read.

- `at`:

  the held coordinate's position in that vector.

- `par`:

  the restricted estimates, stacked.

- `information`:

  a function of no arguments returning the penalized information at the
  restricted point. It is read at most once however often it is called,
  and not at all where nothing calls it, which is the ordinary case: the
  likelihood ratio and the gradient statistic need no matrix.

- `mode_error`:

  a function of no arguments returning
  [`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
  how far above its own mode the refit stopped. It reads the
  information, so it costs one Hessian.

- `labels`:

  the stacked coefficient labels.

- `converged`:

  the inner optimizer's flag, a single logical. Measured, it is
  anti-correlated with how well the point is located – see
  [`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
  which is the reading to prefer.

## Details

The hold is written on the specification and enforced in the solve, so
no formula is rewritten and no term is rebuilt: a design whose bases
were rebuilt would carry different knots and answer for a different
model. See
[`held_positions()`](https://statmodels7.github.io/statmodels7/reference/held_positions.md)
and
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
for the two halves.

## The hyperparameters are held too

The refit runs the alternation ONCE at the hyperparameters the
unrestricted fit arrived at, with no outer search. That makes the
statistics built on it CONDITIONAL on the smoothing the data chose once,
which is the honest reading of a penalized fit and is also what makes an
interval by inversion affordable – re-selecting the smoothing at every
held value would put a whole outer search inside every step of a root
find. A caller who wants the profiled version can refit by hand.

## What is refused

A coefficient under a KINKED penalty – a lasso, a SCAD, an MCP – for two
reasons that point the same way: that block is fitted by a coordinate
descent, and
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
does not read `held_coef` at all, so a hold there would be moved by the
sweep without a word (measured: held at 3, the refit reports 2.7348);
and at the kink the objective has no curvature, so the statistics have
nothing to read and the null distribution of a coefficient a selection
kept is a different problem. An ALIASED coefficient, which has no
estimate to test. And a model carrying a structural term, whose
contribution is a recursion's state rather than \\X\beta\\ and which the
hold does not reach.

A coefficient under a penalty that IS twice differentiable is held like
any other: a smooth's linear column and its rotated coordinates, a
ridge, a random effect. What such a fit means is stated at
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).

## See also

[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md),
which drops the held coordinate from the system.
