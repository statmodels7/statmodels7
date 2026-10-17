# Settle the Sharp Break-Points on the Profile and Hold Them

After the working phase and the restarts, every sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
term is moved to the minimum of the least-squares profile of its working
response
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)),
which is exact only for the location of a gaussian family, held there
([`modelterms7::seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.html)),
and the coefficients of the model are refitted while the positions are
held.

## Usage

``` r
statmod_settle_breakpoints(
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

  The specification, including the holds.

- design:

  The design, carrying the terms in its state.

- blocks, hyper, inner_optimizer, expected, approx, maxit, tol, vb:

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

A list: `res`, the settled fit, `moved`, `TRUE` where a polished
position replaced the one that the iteration reached, and `settled`,
`TRUE` where a term was held.

## Details

The working construction of Fasola, Muggeo and Kuchenhoff stops at a
fixed point of its own iteration. Its coefficients there are those of
the working model, whose block carries a weight frozen at the previous
iterate, and not the coefficients that maximize the likelihood at the
positions reached; and the profile of a discontinuous term is constant
between consecutive observations, so the iteration can stop one interval
away from the minimum. The certificate of such a fit can then read the
mode along the auxiliary column of the working block, which is not a
direction of the model.

Held, a break-point term contributes a linear function of its remaining
coefficients, so the refit is an ordinary fit, and the coordinate that
carried the position is held in the solve and left out of the
information
([`term_held_stack()`](https://statmodels7.github.io/statmodels7/reference/term_held_stack.md)).
`segmented` and `stepmented` report their coefficients in this way, at
the estimated positions and conditional on them.

The profile is weighted least squares of the working response of the
term's own equation, net of the rest of that equation, on the term's own
columns, with the working weights of a scoring step. It is exact for the
location of a gaussian family on the identity link and is a quadratic
model of the objective in any other equation. A polished position is
kept only where the refit at it reaches a better objective than the
refit at the position that the iteration reached. The polish is
repeated, up to `rounds` times, while the refit improves the objective.

A term whose break-point or per-break-point coefficients carry a
development is held and polished like any other term, with the exact
line searches of
[`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)
in the coefficients of the development. In an equation whose working
profile is only a model of the objective, the five deepest local minima
of the profile are also refitted while the positions are held, and so
are the neighboring intervals of the best position; these two searches
skip a developed term, whose positions are one per observation.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which runs before it.
