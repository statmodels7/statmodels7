# What the Fit Certifies About the Point It Reports

Three readings taken AT THE REPORTED POINT and independent of the path
the search took: the outer criterion's gradient, how far the
coefficients sit above the penalized mode, and which hyperparameters
have run to a boundary.

## Usage

``` r
statmod_certificate(fit, tol = 0.01, edge = 8)
```

## Arguments

- fit:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- tol:

  The largest outer gradient a certified point may carry.

- edge:

  The free value beyond which a hyperparameter whose gradient has
  already met `tol` is reported as sitting at a boundary. It decides the
  label alone and never the verdict; see the details.

## Value

A list with `state` (`"converged"`, `"boundary"`, `"not converged"` or
`"unknown"`), `gradient`, `mode_error`, `boundary` and `reason`.

## Details

**Why a certificate rather than the optimizer's flag.** The flag says
whether a search stopped on its own rule, which is a statement about the
search. Measured across shapes, it does not order fits by quality: on
one model the default reported success at a criterion of -1783.47 while
the same data under
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
reached -1664.43 and reported failure. What a reader wants is a property
of the point.

**The state comes from the gradient and the mode error is reported
beside it, not folded into it.** Measured at the reported point over six
shapes, the outer gradient separates by five orders – 4.7e-07, 7.8e-07,
5.8e-05, 7.7e-05 and 3.0e-04 on fits that are right, against 28.8 on one
that is not – while the mode error does not: it reads 1.8e-16 to 6.1e-12
on four of them, 22.8 on the failing one, and 0.114 on a
random-changepoint `seg` whose answer is right to a correlation of
0.9932. A number that does not separate cannot decide a state, and a
certificate that says how far from the mode is worth more than a boolean
that hides it.

`tol` is 1e-2 rather than the geometric middle of the two groups: the
two ways of being wrong are not symmetric, and a certificate that says
NOT CONVERGED at a good point is visible and checkable where one that
certifies a bad point is the failure this exists to remove.

**What it costs** is one outer gradient and one solve, once, at a point
the fit already holds. Nothing is refitted: measured, the criterion
reconstructed from `fit@spec` equals the one the fit reports EXACTLY on
every shape, so the reading is of the fitted model and not of another
one.

**Where there is no outer gradient there are two cases, and they get
different answers.** A model with NO PENALTY – `linpar`, `nl`, `seg`,
`jump`, `jseg` – has no hyperparameter for a gradient to be about, so
the only question left is whether the inner fit reached its mode, and
the mode error answers it: measured over the reference battery it reads
5.2e-11 to 7.9e-05 on fits that are right against 1.215 on a `jump`
fitted to data carrying a slope and a slope change it has no term for. A
model whose only hyperparameters are KINKED – `lasso`, `scad`, `mcp`,
swept along a path because a Laplace approximation at a mode sitting on
the kink has no meaning – gets neither reading and stays `"unknown"`: at
a coefficient the penalty has set to zero the score does not vanish but
lies in the subdifferential, so the mode error is not a statement about
being at a mode. Measured on a lasso, its 4.7e-03 is carried by a
coordinate whose coefficient is exactly 0 and whose score is -0.715.

A form whose criterion has no EXACT gradient
([`outer_gradient_ok`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md))
is also `"unknown"` rather than approximated: 2p refits to difference it
would cost more than the fit.

**The boundary label, and why its threshold needs no derivation.** A
hyperparameter may run to an edge and belong there: on a covariate that
is pure noise the smoothing parameter reaches 9.2e+08, the criterion is
genuinely flat, and calling that fit unconverged would be wrong. A
coordinate is reported as sitting at a boundary when its free value
exceeds `edge` AND its own gradient component has already met `tol`.
Because of that second condition the threshold cannot change the
verdict: a coordinate it moves out of the interior set had already
passed the test, so the maximum that decides the state is unaffected,
and both `"converged"` and `"boundary"` are certified. What `edge`
decides is how the point is described. The default separates the
measured cases with room on both sides: coordinates that ran to an edge
sit at 9.3, 10.5 and 20.6 on the free scale against 0.13, 0.30 and 2.01
for the ones that did not.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`mode_error_limit`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
[`criterion_resolution`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)

## Examples

``` r
dd <- data.frame(x = runif(120))
dd$y <- sin(4 * dd$x) + rnorm(120, 0, 0.3)
statmod_certificate(statmod(y ~ s(x, k = 8),
                            distributions7::gaussian1_distrib(), dd))$state
#> [1] "converged"
```
