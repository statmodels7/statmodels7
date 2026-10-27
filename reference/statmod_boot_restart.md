# Bootstrap Restarting of a Break-Point Fit

Improves a fitted break-point model by bootstrap restarting (Wood 2001),
screening each proposal on the profile of the term's equation before a
refit of the model. The objective of a model with break-points has
spurious local minima, and an ordinary fit converges into the one in
which its starting positions lie; this search looks for a better one.

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
  emethod,
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

  The hyperparameters at which the fit ended, held fixed throughout. The
  restarts search over positions only.

- inner_optimizer:

  How the smooth block is fitted,
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- res:

  The fitted result to improve, as
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  returns it. It is returned unchanged when no proposal improves the
  objective, and also when the design carries no state, or when the
  response is a matrix, is not numeric or has missing values.

- expected, emethod, maxit, tol:

  Passed to
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  for each refit, with the meanings that they have there.

- vb:

  The resolved verbosity, as
  [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  returns it.

- nb:

  The budget: the maximum number of proposals to try, as
  [`seg_boot_total()`](https://statmodels7.github.io/statmodels7/reference/seg_boot_total.md)
  reports it. The first proposal is the deterministic sweep. A budget
  below one returns `res` unchanged.

## Value

`res`, in which `par`, `value`, `converged`, `obj` and the block
histories are replaced when a restart improves the objective. Every
other element is carried over unchanged, including the history of an
outer search and its optimizer, so the result has the same structure
either way.

## The screen on the profile

The non-convexity of a break-point model is in the positions. Given the
positions, the coefficients of a term enter linearly, so the profile at
a fixed configuration of positions is one weighted linear fit of the
working response of the term's equation. A proposal is therefore a
configuration of positions, and two proposals are ranked by their
profiles at the cost of one linear fit each. Only a proposal that ranks
higher on the profile leads to a refit of the whole model, and the refit
is accepted or rejected on the true objective. A proposal that the
profile rejects costs profile evaluations only, with no refit of the
model.

## The three proposal kinds

The kinds are tried in this order.

1.  The deterministic sweep. Each break-point in turn is swept over a
    grid on the profile while the others are held, as
    [`modelterms7::seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.html)
    does. This reaches a feature of the profile from which the fitting
    iteration had moved a break-point away.

2.  The bootstrap sweep. The same descent on the profile of a resample,
    in which the multinomial counts enter as weights. This moves the
    profile's optima in the way that refitting the resample would.

3.  The random sweep. The same descent from positions drawn uniformly
    over the confinement interval.

The two stochastic kinds alternate. Four consecutive proposals that do
not improve the objective, either because the screen rejects them or
because the refit does not decrease the objective, end the loop,
whatever budget is left.

## The profile for a gaussian location and for other equations

The profile is the weighted least-squares fit of the working response of
the term's equation, net of the other contributions to that equation,
with the working weights of a scoring step. For the location parameter
of a gaussian family on the identity link, this is the least-squares
objective of the model. In every other equation, it is a quadratic model
of the objective that ranks the candidates, and each candidate is
accepted or rejected on the true objective.
[`modelterms7::seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
uses the same profile to place starting values.

## Reproducibility

The draws come from the session's generator, so a fit with restarts
repeats under [`set.seed()`](https://rdrr.io/r/base/Random.html). The
refreshable and structural state of the design is saved at the incumbent
and restored whenever a candidate is rejected, so that a rejected
proposal leaves the state unchanged.

## References

Wood, S. N. (2001). Minimizing model fitting objectives that contain
spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
240–244.

## See also

[`seg_boot_total()`](https://statmodels7.github.io/statmodels7/reference/seg_boot_total.md),
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
