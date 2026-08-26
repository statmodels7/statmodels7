# Restarting Around a Fitted Model, Screened on the Exact Profile

Improves a fitted break-point model by bootstrap restarting (Wood 2001),
screening each proposal on the exact profile before paying for a refit.
The objective of a model with break-points has spurious local minima,
and an ordinary fit converges into whichever one its starting positions
sit in; this searches for a better one.

## Usage

``` r
statmod_boot_restart(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  res,
  expected,
  approx,
  maxit,
  tol,
  vb,
  nb
)
```

## Arguments

- spec:

  The
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  being fitted.

- design:

  The assembled design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- blocks:

  The split of the terms into the jointly fitted smooth block and the
  kinked ones, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns it.

- hyper:

  The hyperparameters the fit ended at, held fixed throughout: the
  restarts search over positions, not over hyperparameters.

- inner_optimizer:

  How the smooth block is fitted,
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- res:

  The fitted result to improve, as
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  returns it. Returned unchanged when nothing better is found.

- expected, approx, maxit, tol:

  Passed to
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  for each refit, with the same meanings they have there.

- vb:

  The resolved verbosity, as
  [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  returns it.

- nb:

  The budget: at most how many proposals to try, as
  [`seg_boot_total()`](https://statmodels7.github.io/statmodels7/reference/seg_boot_total.md)
  reports it. A budget of zero returns `res` untouched.

## Value

`res`, with `par`, `value`, `converged`, `obj` and the block histories
replaced when a restart improved the objective. Everything else is
carried over untouched, an outer search's history and its optimizer
among them, so the result is the same shape either way.

## What makes the screen cheap

The non-convexity of a break-point model lives entirely in the
positions. Hold them and everything left is convex, so the exact profile
at a fixed configuration of positions is one linear fit. A proposal is
therefore a configuration of positions and nothing else, and two
proposals are ranked by their profiles at a cost of one linear fit each.
Only a proposal the profile prefers earns a refit of the whole model,
and the refit's answer is accepted or rejected on the true objective.

Measured at \\n = 10^4\\: a proposal that goes nowhere costs about half
a second, against 5 to 15 seconds for the refit it would otherwise have
triggered. The design this replaced refitted every proposal and spent
945 seconds re-verifying an optimum the sweep had already found.

## The three proposal kinds

Tried in this order.

1.  **The deterministic sweep.** Each break-point in turn is swept over
    a grid on the profile with the others held, which is
    [`modelterms7::seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.html).
    This walks straight to a feature the fitting iteration pressed a
    break-point away from.

2.  **The bootstrap sweep.** The same descent on the profile of a
    resample, the multinomial counts entering as weights, which moves
    the profile's optima the way refitting the resample would.

3.  **The random sweep.** The same descent from positions drawn
    uniformly over the confinement interval.

The two stochastic kinds alternate. Four consecutive proposals that fail
the screen end the loop, whatever budget is left.

## The profile is exact for an identity link and a proposal elsewhere

It reads the response net of the other contributions in the term's
equation, on the predictor scale. For a Gaussian response and an
identity link that is the model's own least-squares objective. For
anything else it is an approximation used to rank candidates, and the
true objective decides the acceptance, so a poor ranking costs time and
never correctness.
[`modelterms7::seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
makes the same argument for the same reason.

## Reproducibility

The draws come from the session's generator, so a fit with restarts
repeats under [`set.seed()`](https://rdrr.io/r/base/Random.html). The
refreshable and structural state of the design is snapshotted at the
incumbent and restored whenever a candidate loses, so a rejected
proposal leaves nothing behind.

## References

Wood, S. N. (2001). Minimizing model fitting objectives that contain
spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
240–244.

## See also

[`seg_boot_total()`](https://statmodels7.github.io/statmodels7/reference/seg_boot_total.md),
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
