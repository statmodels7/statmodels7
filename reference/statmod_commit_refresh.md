# Advance the Refresh State

Replaces the terms the design refreshes from by the ones they refresh to
at the given coefficients, so that whatever a term carries as a state of
its iteration moves on one step.

## Usage

``` r
statmod_commit_refresh(spec, coef, design, which = "all")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors, one per distribution parameter.

- design:

  The design, whose refresh state is what this advances.

- which:

  Which entries to commit: `"all"` (the default), `"jacobian"` or
  `"frozen"`. The alternation's pass level and
  [`statmod_fitted_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_fitted_spec.md)
  both pass `"jacobian"`, the frozen terms having been committed already
  by their own phase.

## Value

The coefficient list, invisibly, with each committed term's stretch
replaced by the coefficients that term stored. Identical to `coef` when
nothing was committed or when no term relabeled.

## What moves

The rescaling factor of a discontinuous break-point term, and the
direction it last traveled in, which modelterms7 halves on a reversal.
Advancing that once per objective evaluation would anneal at the speed
of the line search instead of the speed of the fit; never advancing it
would solve a permanently smoothed problem, whose fixed point is not the
model's. Once per sweep is what this call is for.

## Why a frozen block is committed elsewhere

Where a term's block is a Jacobian, committing does not move the
objective at the same coefficients: the break-point is read off the
coefficients and the rescaling reaches only the columns.

Where the block is a frozen working linearization, that is false twice
over, which is why those terms are committed by
[`fit_working()`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
and skipped here at the default. A
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
reads its position from a quadratic that is incremental in the position
already committed, so a second commit at the same coefficients takes a
further hidden step, measured at up to 0.71 per observation on the
contribution. And a refresh may relabel crossed break-point lineages,
after which the caller's own copy of the coefficients names them in the
old order.

The relabeling is why the committed coefficients are returned. A caller
continues from what the terms stored, not from what it passed in.

## See also

[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the refresh this advances,
[`fit_working()`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
for the phase that commits the frozen terms,
[`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md)
for the verdict.
