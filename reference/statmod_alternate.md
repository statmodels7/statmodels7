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

  The block split.

- hyper:

  The hyperparameters.

- inner_optimizer:

  How the smooth block is fitted.

- beta:

  The starting coefficients, stacked.

- expected:

  Whether the information is the expected one.

- approx:

  The approximation for the expected information.

- maxit, tol:

  The budget and the tolerance.

- vb:

  The resolved verbosity.

- working_budget:

  How many working fits
  [`fit_working`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
  may take. The bootstrap excursions of
  [`statmod_boot_restart`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
  pass a short one: an excursion needs to travel, not to converge.

- hold_refresh:

  `TRUE` holds every FROZEN break-point block at its committed
  positions: the fit is then an ordinary smooth fit, no read-off runs
  and no schedule advances. The outer machinery passes it at every
  criterion evaluation, path point and fold, for two reasons measured
  together: the working phase inside each of dozens of evaluations
  multiplied a fit's cost by twenty, and a break-point moving between
  evaluations makes the criterion path-dependent while its cycling flags
  read as unavailable points to the search. The positions are refined
  ONCE, by the full alternation
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  runs at the chosen hyperparameters before the restarts.

## Value

A list with `par`, `value`, `converged`, `obj`, `hist_blocks` and
`hist_inner`.

## Details

It is a function of its own because the outer search calls it once per
hyperparameter it tries, warm-started from the previous coefficients.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
