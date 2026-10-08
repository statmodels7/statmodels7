# Settle the Sharp Break-Points on the Profile and Hold Them

After the working phase and the restarts, every sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
term is moved to the minimum of its exact least-squares profile
([`modelterms7::seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.html)),
held there
([`modelterms7::seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.html)),
and the coefficients of the model are refitted with the positions held.

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

A list: `res`, the settled fit, `moved`, `TRUE` where a polished
position replaced the one the iteration reached, and `settled`, `TRUE`
where a term was held.

## Details

The working construction of Fasola, Muggeo and Kuchenhoff stops at a
fixed point of its own iteration. Its coefficients there are those of
the working model, whose block carries a weight frozen at the previous
iterate, and not the coefficients that maximize the likelihood at the
positions reached; and the profile of a discontinuous term is constant
between consecutive observations, so the iteration can stop one interval
away from the minimum. Measured on 16 samples of 200 observations (a
jseg and a jump, 8 seeds each), the fit reached the interval of the
profile's minimum in 11, its coefficients sat up to 0.04 log-likelihood
units below least squares at the position reached, and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reported `not converged` in 12, reading the mode along the working
block's auxiliary column, which is not a direction of the model.

Held, a break-point term contributes a linear function of its remaining
coefficients, so the refit is an ordinary fit, and the coordinate that
carried the position is held in the solve and left out of the
information
([`term_held_stack()`](https://statmodels7.github.io/statmodels7/reference/term_held_stack.md)).
This is how `segmented` and `stepmented` report their coefficients: at
the estimated positions, conditional on them.

The profile is weighted least squares of the working response of the
term's own equation, net of the rest of that equation, on the term's own
columns, with the working weights of a scoring step: exact for a
gaussian mean and the quadratic model of the objective for any other
equation. A polished position is kept only where the refit at it reaches
a better objective than the refit at the position the iteration reached.
A term whose per-break-point coefficients carry a development is held
without polishing, its positions being one per observation.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which runs before it.
