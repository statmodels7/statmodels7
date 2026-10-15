# Fit the Smooth Block Around a Frozen Working Linearization

The iteration of Fasola, Muggeo and Kuchenhoff (2018) for a term whose
block is a working linearization with a frozen weight
([`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)).
The smooth block is fitted exactly at the committed block, the
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

  How many working fits at most.

- stall_limit:

  How many consecutive working fits with every unsettled break-point at
  its scaling floor end the phase.

## Value

A list with the elements of
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
where `iterations` is the total over the working fits and `history` is
`NULL`, and with `fasola`, the number of working fits taken, and
`stalled`, `TRUE` where the phase was ended by `stall_limit`. In the
iteration proper the list also carries a `note` element, which is
`NULL`.

## Details

The sequencing is the whole of the difference from
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
and it is the sequencing of the segmented package. The fixed-point
iteration these constructions belong to is not a descent method on the
objective of the model: its early steps under a large scaling factor
move uphill on purpose, which is how it leaves a spurious optimum.
Embedding the read-off inside the objective of the inner optimizer would
put a sufficient-decrease line search in its way and stall it. During
the working fit the frozen blocks contribute \\X\beta\\ and nothing else
(`st$working`), which makes the inner fit the plain penalized working
fit of the papers; the commit then advances the read-off, the scaling
schedule and any relabeling of crossed break-point lineages, once per
working fit.

Any inner method serves: the working fit goes through
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which takes
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
or any optimizers7 optimizer, and the read-off never moves inside the
objective of the inner method. The methods differ in cost. Each working
fit is solved afresh at a frozen block, so a method that carries exact
curvature closes it in a step or two, while a quasi-Newton method
rebuilds its own curvature from nothing every time.

The exit is at a fixed point of the iteration or in the cycle it settles
into, judged on the working objective. The iteration stops when the
read-off has settled
([`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html))
with the objective stalled, when the objective has stalled three times
in a row, or when the objective equals its value two iterations back
twice (the period-two cycle of the break-point that Muggeo documents,
which a consecutive-change rule never sees). Only at a fixed point does
the working objective coincide with the model's, so no best-so-far
iterate is kept: mid-travel, the committed contribution of a good
working value can sit orders of magnitude off the data. Running out of
the budget reports `FALSE`.

A phase that cannot settle is ended early and reports `FALSE` as well.
This happens where every break-point that has not settled sits at the
floor of its scaling factor
([`modelterms7::term_stalled()`](https://statmodels7.github.io/modelterms7/reference/term_stalled.html))
for `stall_limit` consecutive working fits: the factor can shrink no
further, and nothing is left to bring the break-point to rest. The
profile objective of a discontinuous term is constant between
consecutive observations, so such a break-point passes from one
observation to the next for as long as the budget lasts. The default
`stall_limit` of 20 ends such a phase, and the stall and cycle rules are
read first at every working fit.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
