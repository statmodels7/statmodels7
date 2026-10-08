# Settle the Changes of Slope on the Profile

After the working phase and the restarts, every sharp
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term is moved to the minimum of its exact least-squares profile over
every interval of the covariate
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)),
and the model is refitted from there with the positions free.

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
  approx,
  maxit,
  tol,
  vb,
  rounds = 5L
)
```

## Arguments

- spec:

  The specification, holds included.

- design:

  The design, carrying the terms in its state.

- blocks, hyper, inner_optimizer, expected, approx, maxit, tol, vb:

  As
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md).

- res:

  The fit to settle, as
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  returns it.

## Value

A list: `res`, the fit, and `moved`, `TRUE` where a polished position
replaced the one the iteration reached.

## Details

The profile of a change of slope is smooth inside each interval between
consecutive observed values and has a kink at each of them, so the
iteration of Muggeo (2003) can stop on a kink that is not a minimum:
measured on
[`segmented::globTempAnom`](https://rdrr.io/pkg/segmented/man/globTempAnom.html)
with four break-points, it stopped with two positions on observed years
and the inner fit 0.106 log-likelihood units above its mode, and the
REML criterion could not be evaluated at its start. Inside an interval
the profile has one stationary point, so its minimum over every interval
is exact
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)).
The profile is read on the working response of the term's equation, and
a polished position is kept only where the refit from it reaches a
better objective. In an equation whose working profile only approximates
the objective, the five deepest local minima of the profile are refitted
as well.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22, 3055–3071.
