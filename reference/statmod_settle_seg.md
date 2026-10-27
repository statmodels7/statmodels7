# Settle the Changes of Slope on the Profile

After the working phase and the restarts, every
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term that is not smoothed, not developed and not marginal is moved to
the minimum of the least-squares profile of its working response over
every interval of the covariate
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)),
and the model is refitted from there while the positions are free.

## Usage

``` r
statmod_settle_seg(
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
  rounds = 5L
)
```

## Arguments

- spec:

  The specification, including the holds.

- design:

  The design, carrying the terms in its state.

- blocks, hyper, inner_optimizer, expected, emethod, maxit, tol, vb:

  As
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md).

- res:

  The fit to settle, as
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  returns it.

- rounds:

  The maximum number of times that the polished positions are read again
  after a refit, 5 by default.

## Value

A list: `res`, the fit, and `moved`, `TRUE` where a polished position
replaced the one that the iteration reached.

## Details

The profile of a change of slope is smooth inside each interval between
consecutive observed values and has a kink at each of them, so the
iteration of Muggeo (2003) can stop on a kink that is not a minimum.
Inside an interval the profile has one stationary point, so its minimum
over every interval is exact
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)).
The profile is read on the working response of the term's equation, and
a polished position is kept only where the refit from it reaches a
better objective. In an equation whose working profile only approximates
the objective, the five deepest local minima of the profile are refitted
as well.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22, 3055–3071.
