# Fit the Smooth Block Around a Frozen Working Linearization

The iteration of fasola2018 for a term whose block is a working
linearization with a frozen weight –
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html):
the smooth block is fitted exactly at the committed block, the
break-points are read off the fitted coefficients and committed, and the
two alternate until the read-off settles or the working objective stops
moving.

## Usage

``` r
fit_working(
  obj,
  beta,
  idx,
  spec,
  design,
  hyper,
  method,
  vb,
  tol,
  budget = 500L,
  stall_limit = 20L
)
```

## Arguments

- obj:

  The objective.

- beta:

  The current stacked coefficients.

- idx:

  The smooth block's indices.

- spec:

  The specification.

- design:

  The design.

- hyper:

  The hyperparameters.

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizer.

- vb:

  The resolved verbosity.

- tol:

  The alternation's tolerance, read for the objective-stall rule.

- budget:

  How many working fits at most. The default covers the measured runs
  (69 to 165 iterations on three break-points) with room.

- stall_limit:

  How many consecutive working fits with every unsettled break-point at
  its scaling floor end the phase.

## Value

As
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
plus `fasola`, the number of working fits taken, and `stalled`, `TRUE`
where the phase was ended by `stall_limit`.

## Details

The sequencing is the whole of the difference from
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
and it is what `segmented` does. The fixed-point iteration these
constructions belong to is not a descent method on the model's objective
– its early steps under a large scaling factor move uphill on purpose,
which is how it leaves a spurious optimum – so embedding the read-off
inside the inner optimizer's objective put a sufficient-decrease line
search in its way and stalled it: measured on a three-break-point jseg,
the embedded route ended at an rss worse than the mean-only fit from the
TRUE break-points, while this iteration recovers them from the same
start. During the working fit the frozen blocks contribute \\X\beta\\
and nothing else (`st$working`), which makes the inner fit the plain
penalized working fit of the papers; the commit then advances the
read-off, the scaling schedule and any relabeling of crossed break-point
lineages, once per working fit.

Any inner method serves: the working fit goes through
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which takes
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
or any optimizers7 optimizer, and the read-off never moves inside
anyone's objective. What differs is the price. Each working fit is
solved afresh at a frozen block, so a method carrying exact curvature
closes it in a step or two while a quasi-Newton method rebuilds its own
from nothing every time: measured on three break-points at \\n =
10000\\, the same answer to the digit costs 6.9 s under
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
9.2 s under `newton()` and 140 s under `lbfgs()`.

The exit is at a fixed point of the iteration or in the cycle it settles
into, judged on the working objective: the read-off settled
([`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html))
with the objective stalled, the objective stalled three times in a row,
or the objective equal to two iterations back twice – the period-two
cycle of the break-point Muggeo documents, which a consecutive-change
rule never sees. Only at a fixed point does the working objective
coincide with the model's, which is why no best-so-far iterate is kept:
mid-travel the committed contribution of a good working value can sit
orders of magnitude off the data. Running out of the budget reports
`FALSE`.

A phase that cannot settle is ended early and reports `FALSE` as well:
where every break-point that has not settled sits at the floor of its
scaling factor
([`modelterms7::term_stalled()`](https://statmodels7.github.io/modelterms7/reference/term_stalled.html))
for `stall_limit` consecutive working fits, the factor can shrink no
further and nothing is left to bring the break-point to rest. The
profile objective of a discontinuous term is constant between
consecutive observations, so such a break-point passes from one
observation to the next for as long as the budget lasts: measured on
`jseg(x, psi = 5, n_boot = 0)` at \\n = 400\\, the factor reached its
floor of 1.8e-8 in the second pass and the break-point then wandered
through a hundred passes of 500 working fits, 245 s, without converging.
Four of those passes did end on the stall or cycle rule, after 58, 116,
349 and 349 working fits at the floor, which is a wandering break-point
landing on a step small enough by chance and not a settled one. In the
phases that settle the longest run at the floor measured is 7 working
fits (three break-points at \\n = 10000\\); the default of 20 sits
between the two, and the stall and cycle rules are read first at every
working fit. With it the fit above ends in 1.7 s, and with the restarts
left on it reaches the same optimum as the grid start in 1.3 s where it
took 250 s.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
