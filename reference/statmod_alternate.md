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
  emethod,
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

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

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
  positions: the fit is then an ordinary smooth fit, and the schedule is
  neither read off nor advanced. The outer machinery passes it at every
  criterion evaluation, path point and fold. A working phase inside each
  evaluation would multiply the cost of the fit, and a break-point that
  moves between evaluations makes the criterion depend on the path of
  the search. The positions are refined once, by the full alternation
  that
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  runs at the chosen hyperparameters before the restarts.

## Value

A list of nine:

- `par`:

  the stacked coefficients reached.

- `value`:

  the penalized objective there, not averaged.

- `converged`:

  a single logical. Where a term recomputes its own block, the verdict
  is the objective's together with that term's own flags, and the flag
  of the smooth block is not part of it. Beside a kinked block or a
  structural term it also requires
  [`alternation_settled()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md)
  at the point. Where the fit does not alternate between blocks, it is
  the flag of the inner method.

- `obj`:

  the objective object that the last pass used.

- `note`:

  the `note` element of the last smooth fit, as
  [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  returns it.

- `aliased`:

  the sorted stacked positions that the pivot left out.

- `hist_blocks`:

  a list of one-row data frames, one for each block fitted in each pass,
  with `pass`, `block`, `objective`, `change`, `iterations` and
  `converged`.

- `hist_inner`:

  a list of data frames, one for each pass in which the smooth fit
  recorded steps, holding the inner history with the columns `pass` and
  `block` added.

- `spec`:

  the specification after the last pass, with the SCAD and MCP
  curvatures that the fit used.

## Details

The outer search calls this function once for each hyperparameter value
that it tries, warm-started from the previous coefficients, so a whole
fit at one point of the search is one call.

One pass fits the smooth block (beside a structural term of the filter
shape, the joint system of the coefficients and the term's parameters),
then the structural block of a term of the likelihood shape, then each
kinked block while the other blocks are held, and then commits the
schedule of any refreshable term. The passes stop when the relative
change of the objective falls below `tol`. The verdict is the inner
fit's flag where the fit does not alternate between blocks. Otherwise it
also requires
[`alternation_settled()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md)
(the Newton decrement of the smooth part and the KKT decrement of each
coefficient that a kinked penalty holds at zero, both within
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)),
and
[`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md)
where the model has refreshable terms.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the fit that this performs,
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
and
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
for the two halves of a pass, and
[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which calls this repeatedly.
