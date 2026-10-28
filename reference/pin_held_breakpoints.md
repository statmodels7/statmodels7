# Keep the Free Break-Points Away From a Held One

Prepares the break-point terms of a restricted fit in which a
coefficient of a continuous
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term with two or more break-points is held. The free break-points are
confined to the side of the held position on which they start, with at
least `min_obs` observations between them and it
([`modelterms7::seg_pin()`](https://statmodels7.github.io/modelterms7/reference/seg_pin.html)),
and the starting values of their positions are brought inside the
interval. The position of a free break-point is kept inside its interval
by an infinite objective outside it
([`pin_infeasible()`](https://statmodels7.github.io/statmodels7/reference/pin_infeasible.md)).

## Usage

``` r
pin_held_breakpoints(
  spec,
  design,
  beta,
  from = beta,
  min_obs = 2L,
  variant = 1L
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  with the coefficient held.

- design:

  Its design, whose state holds the terms that a fit refreshes. The
  terms in the state are replaced.

- beta:

  The starting coefficients, stacked, with the held value in place.

- from:

  The coefficients of the unrestricted fit, stacked, from which the
  starting pairs of position and change of slope are taken.

- min_obs:

  The smallest number of observations between a free and a held
  break-point.

- variant:

  Which unrestricted break-point the held position stands for, counted
  from the nearest to it: the held position takes the change of slope of
  that break-point and the free ones keep the other pairs.

## Value

A list with `beta`, the starting coefficients with the free positions
inside their intervals, `present`, a function of the coefficient list of
the restricted fit that returns it with the break-points in order, and
`variants`, the number of break-points that the held position can stand
for: the number of break-points of the term when one position is held,
and one otherwise. A model without such a term is left unchanged.

## Details

Two break-points at one position make two columns of the block equal. A
free break-point that reaches the held one then has coefficients of the
changes of slope that grow without bound in opposite directions, and the
relabeling of crossed break-points of
[`modelterms7::term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html)
writes the held value onto the coefficients of the free one. A
restricted fit of this kind was measured at 4006 scoring iterations,
without converging, for a statistic that was not the likelihood ratio.

The coefficients of such a term are returned in the order of the
positions of the break-points, by `present()`, and the coefficient
labeled as held carries `[held]` after its name.

## See also

[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
which calls it.
