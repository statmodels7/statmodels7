# The Alternation Between the Smooth Block and the Rest

Fits the terms whose penalties are twice differentiable in one system
and each remaining block by a method of its own, alternating until the
objective stops moving.

## Usage

``` r
statmod_alternate(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  beta,
  expected,
  approx,
  maxit,
  tol,
  vb,
  working_budget = 500L,
  hold_refresh = FALSE
)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- blocks:

  The split of the terms into `smooth` and `sparse`, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns it.

- hyper:

  The hyperparameters, held fixed for the whole call.

- inner_optimizer:

  How the smooth block is fitted:
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- beta:

  The starting coefficients, stacked into one numeric vector.

- expected:

  `TRUE` for the expected information, `FALSE` for the observed one.

- approx:

  How the expected information is approximated for a family with no
  closed form.

- maxit, tol:

  The budget in passes and the relative tolerance on the objective, as
  [`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
  read them off `inner_optimizer`.

- vb:

  The resolved verbosity, as
  [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  returns it.

- working_budget:

  How many working fits
  [`fit_working()`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
  may take. The bootstrap excursions of
  [`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
  pass a short one: an excursion needs to travel, not to converge.

- hold_refresh:

  `TRUE` holds every frozen break-point block at its committed
  positions: the fit is then an ordinary smooth fit, no read-off runs
  and no schedule advances. The outer machinery passes it at every
  criterion evaluation, path point and fold, for two reasons measured
  together: the working phase inside each of dozens of evaluations
  multiplied a fit's cost by twenty, and a break-point moving between
  evaluations makes the criterion path-dependent while its cycling flags
  read as unavailable points to the search. The positions are refined
  once, by the full alternation
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  runs at the chosen hyperparameters before the restarts.

## Value

A list of six:

- `par`:

  the stacked coefficients reached.

- `value`:

  the penalized objective there, unaveraged.

- `converged`:

  a single logical: the objective settled, every block settled, and
  every refreshable term reported it had.

- `obj`:

  the objective object the last pass used.

- `hist_blocks`:

  a data frame, one row per pass.

- `hist_inner`:

  a data frame, one row per inner iteration.

## Details

This is a function of its own because the outer search calls it once per
hyperparameter it tries, warm-started from the previous coefficients. A
whole fit at one point of the search is one call of this.

One pass does the smooth block, then each kinked block with everything
else held, then commits any refreshable term's schedule. The passes
repeat until the objective's relative change falls below `tol`, every
block reports it has settled, and
[`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md)
agrees.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the fit this performs,
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
and
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
for the two halves of a pass,
[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which calls this repeatedly.
