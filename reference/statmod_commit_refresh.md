# Advance the Refresh State

Replaces the terms the design refreshes from by the ones they refresh to
at the given coefficients, so that whatever a term carries as a state of
its iteration moves on one step.

## Usage

``` r
statmod_commit_refresh(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

## Value

`TRUE` when there was something to commit, invisibly.

## Details

What moves is the rescaling factor of a discontinuous break-point term
and the direction it last travelled in, which modelterms7 halves on a
reversal: a schedule that advanced once per objective evaluation would
anneal at the speed of the line search rather than at the speed of the
fit, and one that never advanced would solve a permanently smoothed
problem, whose fixed point is not the model's.

The value a term reports is unchanged by the schedule – a break-point is
read off the coefficients and the rescaling reaches only the columns –
so committing does not move the objective at the same coefficients.

## See also

[`statmod_design_at`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
