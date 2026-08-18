# Restarting Around a Fitted Model, Screened on the Exact Profile

Bootstrap restarting (Wood 2001) with the observation that makes it
cheap: the non-convexity of a break-point model lives entirely in the
positions, everything else being convex once they are held, so a
proposal is a CONFIGURATION OF POSITIONS and two configurations are
compared on the exact profile – least squares at fixed positions, one
linear fit each – rather than by refitting the model. Only a proposal
the profile prefers earns a refit, and the refit's answer is accepted on
the true objective.

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

  The specification.

- design:

  The design.

- blocks:

  The block split.

- hyper:

  The hyperparameters the fit ended at.

- inner_optimizer:

  How the smooth block is fitted.

- res:

  The fitted result the restarts try to improve.

- expected, approx, maxit, tol:

  As in
  [`statmod_alternate`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md).

- vb:

  The resolved verbosity.

- nb:

  At most how many proposals.

## Value

`res`, with `par`, `value`, `converged`, `obj` and the block histories
replaced when a restart improved the fit; any other field – an outer
search's history, its optimizer – is kept.

## Details

Three proposal kinds, in order. The first is deterministic: each
break-point swept over the profile with the others held
([`seg_polish`](https://statmodels7.github.io/modelterms7/reference/seg_polish.html)),
which walks straight to a feature the iteration pressed a break-point
away from. The stochastic ones alternate a BOOTSTRAP sweep – the same
descent on the profile of a resample, the multinomial counts entering as
weights, which moves the profile's optima the way refitting the resample
would – with a sweep from RANDOM positions drawn over the confinement
interval. The screen compares every proposal's unweighted profile
against the incumbent's, so a proposal that lands back on the incumbent
costs a few grid-many linear fits and no refit: measured at n = 10000, a
dry proposal fell from a whole refit (5 to 15 s) to about half a second,
and the earlier design that refitted every proposal had spent 945 s
re-verifying an optimum the sweep had already found. Four consecutive
dry proposals end the loop.

The profile reads the response net of the other contributions in the
term's equation, on the predictor scale: exact for an identity link, a
proposal elsewhere – the argument
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
already makes. The refreshable and structural state of the design is
snapshotted at the incumbent and restored whenever a candidate loses,
and the draws come from the session's generator, so a fit with restarts
is reproducible under
[`set.seed()`](https://rdrr.io/r/base/Random.html).

## References

Wood, S. N. (2001). Minimizing model fitting objectives that contain
spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
240–244.

## See also

[`seg_boot_total`](https://statmodels7.github.io/statmodels7/reference/seg_boot_total.md),
[`statmod_alternate`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
