# Refit With One Coefficient Held at a Value

Maximizes the penalized likelihood over every coefficient but one, which
is held at `value`. It returns the restricted maximum, the coefficients
that attain it and the score there.

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

  The value at which to hold it, a single number.

## Value

A list:

- `coefficients`:

  the restricted estimates, a named list.

- `loglik`:

  the log-likelihood there, unpenalized.

- `objective`:

  the penalized objective there, the negative log-likelihood plus the
  penalty, which the fit minimizes.

- `score`:

  the penalized score at the restricted point, a numeric vector over the
  stacked coefficients. Its free entries are zero to the tolerance; the
  statistics read the held one.

- `at`:

  the held coordinate's position in that vector.

- `par`:

  the restricted estimates, stacked.

- `information`:

  a function of no arguments returning the penalized information at the
  restricted point. It is read at most once however often it is called,
  and it is not read if the function is never called, which is the
  ordinary case: the likelihood ratio and the gradient statistic need no
  matrix.

- `mode_error`:

  a function of no arguments returning
  [`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
  how far above its own mode the refit stopped. It reads the
  information, so it costs one Hessian.

- `labels`:

  the stacked coefficient labels.

- `converged`:

  the inner optimizer's flag, a single logical. It records whether a
  stopping rule fired and does not measure the distance to the mode,
  which
  [`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md)
  measures.

## Details

The hold is written on the specification and enforced in the solve, so
the refit does not rewrite the formula and does not rebuild any term: a
design whose bases were rebuilt would carry different knots and would
correspond to a different model. See
[`held_positions()`](https://statmodels7.github.io/statmodels7/reference/held_positions.md)
and
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
for the two halves.

## The hyperparameters are held too

The refit runs the alternation once at the hyperparameters of the
unrestricted fit, with no outer search. The statistics built on it are
therefore conditional on the smoothing that the data chose once. That is
the reading that a penalized fit supports, because the penalty is a
log-prior at the hyperparameters in force, and it keeps an interval by
inversion affordable: re-selecting the smoothing at every held value
would put a whole outer search inside every step of a root find. A
profiled version requires refitting by hand.

## Coefficients that cannot be held

A coefficient under a kinked penalty (a lasso, a SCAD, an MCP) cannot be
held, for two reasons. The block is fitted by a coordinate descent,
which updates every coordinate of the block and cannot leave one at a
held value, so that
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
signals an error when a held coordinate lies in its block. And at the
kink the objective has no curvature for the statistics to read, while
the null distribution of a coefficient that a selection kept is a
different problem. An aliased coefficient has no estimate to test. A
model that carries a structural term cannot hold a coefficient either,
because the contribution of the term is the state of a recursion and not
\\X\beta\\, and the hold does not reach it.

A coefficient under a penalty that is twice differentiable is held like
any other: a smooth's linear column and its rotated coordinates, a
ridge, a random effect.
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md)
states the meaning of such a fit.

## See also

[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md),
which drops the held coordinate from the system.
