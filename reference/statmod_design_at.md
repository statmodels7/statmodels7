# The Design at Given Coefficients

The design with every refreshable term's block recomputed at the
coefficients it currently holds, and the difference between its
contribution and its linearization carried as an adjustment to the
predictor.

## Usage

``` r
statmod_design_at(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design, as
  [`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A design.

## Details

A design with no refreshable term is returned unchanged, so a model of
ordinary terms pays nothing and reaches exactly the same arithmetic as
before.

The refresh is CHAINED from the term the state holds rather than from
the specification, because the rescaling factor of a discontinuous
break-point term is a state of the iteration and not a function of the
coefficients: it halves when the break-point reverses direction, which
is a fact about the path and not about the point. The state advances
only when
[`statmod_commit_refresh`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
is called, so the trial points of a line search all see the same
schedule and the schedule advances once per sweep.

The result is memoized on the coefficients, since the objective, its
gradient and its curvature are asked for at the same point in turn.

## See also

[`statmod_commit_refresh`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
