# The Design at Given Coefficients

Recomputes every refreshable term's design block at the coefficients
currently in hand, and carries the difference between what the term
contributes and what its linearization contributes as a per-observation
adjustment to the predictor. With it, every cross-product in the
objective reads the block as an ordinary design while the predictor
stays exact:

\$\$\mathrm{adj} = \mathrm{term\\value}(\beta) - X(\beta)\\\beta\$\$

For a term whose block is a Jacobian,
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
and
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
the adjustment is non-zero and the resulting scoring step is the
Gauss-Newton one. For
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
the columns satisfy \\X\beta = \mathrm{value}\\ exactly, so the
adjustment is zero and the same step is Fasola's fixed-point iteration.

## Usage

``` r
statmod_design_at(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors, one per distribution parameter,
  each as long as its equation's design is wide.

- design:

  The design to refresh, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A design of the same shape as `design`, with the refreshable terms'
column blocks replaced and each equation's `adj` set to the
per-observation adjustment above. `design` itself when nothing
refreshes.

## A model with no refreshable term pays nothing

[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md)
is empty there, the design is returned as it arrived, and the arithmetic
downstream is untouched.

## Chained from the term, not from the specification

The refresh reads the term the design state currently holds, never the
one the specification was built with. A discontinuous break-point term
carries a rescaling factor that is a state of its iteration: modelterms7
halves it whenever the break-point reverses direction, so it records the
path taken and cannot be recovered from the point reached. Refreshing
from the specification each time would reset that factor to its starting
value and solve a permanently smoothed problem, whose fixed point is not
the model's.

The state advances only when
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
is called, so every trial point of a line search sees one schedule and
the schedule advances once per sweep.

## Memoization

The result is cached on the coefficients, because the objective, its
gradient and its curvature are all asked for at the same point in turn
and each would otherwise rebuild the same blocks.

## See also

[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
to advance the state afterwards,
[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md)
for the terms this walks,
[`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md)
for the verdict.
