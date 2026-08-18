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
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- which:

  Which refresh entries to commit: `"all"`, `"jacobian"` (the default at
  the alternation's pass level, where the frozen ones are already
  committed by their own phase) or `"frozen"`.

## Value

The coefficient list, with each committed term's stretch replaced by the
coefficients the term stored, invisibly.

## Details

What moves is the rescaling factor of a discontinuous break-point term
and the direction it last travelled in, which modelterms7 halves on a
reversal: a schedule that advanced once per objective evaluation would
anneal at the speed of the line search rather than at the speed of the
fit, and one that never advanced would solve a permanently smoothed
problem, whose fixed point is not the model's.

For a term whose block is a Jacobian the value it reports is unchanged
by the schedule – a break-point is read off the coefficients and the
rescaling reaches only the columns – so committing does not move the
objective at the same coefficients. For a FROZEN working block that
sentence is false in two ways, which is why those terms are committed by
[`fit_working`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
and skipped here: a jseg's quadratic read-off is incremental in the
committed position, so a second commit at the same coefficients takes a
second step, and a refresh may relabel crossed break-point lineages,
after which the caller's coefficients are stale. The relabeling is why
the COMMITTED coefficients are returned: a caller continues from what
the terms stored, not from what it passed in.

## See also

[`statmod_design_at`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md),
[`fit_working`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
