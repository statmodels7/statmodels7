# What the Fit Certifies About the Point It Reports

Three readings taken at the reported point and independent of the path
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
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- tol:

  The largest rise, in the criterion's own units, that a certified point
  may still have available to it. ⚠️ Until 0.127.0 this was a threshold
  on the outer gradient. The default has not moved and its meaning has:
  a caller who set one should read
  [`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md).

- edge:

  The free value beyond which a hyperparameter that has already met
  `tol` on its own is reported as sitting at a boundary. It decides the
  label alone and never the verdict; see the details.

## Value

A list with `state` (`"converged"`, `"boundary"`, `"not converged"` or
`"unknown"`), `decrement`, `gradient`, `mode_error`, `curvature`,
`boundary`, `boundary_key` and `reason`. `decrement` is what the verdict
is made on and `gradient` is reported beside it; `curvature` says
whether that decrement was read against an analytic Hessian or against
one differenced from the exact gradient. `boundary_key` names the same
coordinates as `boundary` does, in the key
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
labels its rows by, which is what
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reads to leave a standard error off a coordinate pinned there.

## Details

**Why a certificate and not the optimizer's flag.** The flag says
whether a search stopped on its own rule, which is a statement about the
search. Measured across shapes, it does not order fits by quality: on
one model the default reported success at a criterion of -1783.47 while
the same data under
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
reached -1664.43 and reported failure. What a reader wants is a property
of the point.

**The state comes from the rise the criterion would still buy, and the
mode error is reported beside it rather than folded into it.** That rise
is the Newton decrement \\\tfrac{1}{2}g^\top A^{-1}g\\ of
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md),
in the criterion's own units, and `tol` is a threshold on it. The mode
error does not decide a state and could not: measured over six shapes it
reads 1.8e-16 to 6.1e-12 on four fits that are right, 22.8 on one that
is not, and 0.114 on a random-changepoint `seg` whose answer is right to
a correlation of 0.9932. A number that does not separate cannot decide a
state, and a certificate that says how far from the mode is worth more
than a boolean that hides it.

**A threshold on the GRADIENT was what this read until 0.127.0, and it
is not scale-free.** The gradient of a criterion summed over \\n\\
observations and \\p\\ penalized coefficients carries both, so one
number cannot serve every shape: measured over 1350 fits against an
independently located optimum, `max|g| > 1e-2` flags **131 of 663 fits
that are within 1e-3 of their own optimum** while finding all 552 that
are more than 1e-2 short of it. The decrement at the same cut flags **0
of the 663** and finds the same 552. See
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
for the calibration behind that and for the two scaled gradients that
were measured and rejected.

`tol` is 1e-2 rather than the middle of the two groups: the two ways of
being wrong are not symmetric, and a certificate that says not converged
at a good point is visible and checkable where one that certifies a bad
point is the failure this exists to remove.

**What it costs** is one outer gradient, one curvature and one solve, at
a point the fit already holds; the criterion reconstructed from
`fit@spec` equals the one the fit reports exactly on every shape, so the
reading is of the fitted model and of no other. Where the form carries
an analytic outer Hessian nothing is refitted. Where it does not – a
criterion asked for on the expected information, or a separable penalty
– the curvature comes from one central difference of the **exact**
gradient, which is \\4n_h\\ refits and was measured at 0.05 to 0.13
seconds, 5 to 37 per cent of the fit itself.
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md)
says which route was taken, the result reports it in `curvature`, and
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
says so on the certificate's own line where it was differenced – a
verdict resting on a curvature the package estimated rather than
computed should not have to be asked for.

**Where there is no outer gradient there are two cases, and they get
different answers.** A model with **no penalty**, which covers `linpar`,
`nl`, `seg`, `jump` and `jseg`, has no hyperparameter for a gradient to
be about, so the only question left is whether the inner fit reached its
mode, and the mode error answers it: measured over the reference battery
it reads 5.2e-11 to 7.9e-05 on fits that are right against 1.215 on a
`jump` fitted to data carrying a slope and a slope change it has no term
for. A model whose only hyperparameters are **kinked**, `lasso`, `scad`,
`mcp`, swept along a path because a Laplace a Laplace approximation at a
mode sitting on the kink having no meaning, gets neither reading and
stays `"unknown"`: at a coefficient the penalty has set to zero the
score does not vanish but lies in the subdifferential, so the mode error
is not a statement about being at a mode. Measured on a lasso, its
4.7e-03 is carried by a coordinate whose coefficient is exactly 0 and
whose score is -0.715.

A form whose criterion has no exact GRADIENT
([`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
at order 1) is also `"unknown"` and never approximated: 2p refits to
difference it would cost more than the fit. That is not in tension with
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md)
differencing a HESSIAN where the analytic one is missing – what it
differences there is the exact gradient, so the reading still rests on a
derivative the package computes rather than on one it estimates twice
over.

**The boundary label, and why its threshold needs no derivation.** A
hyperparameter may run to an edge and belong there: on a covariate that
is pure noise the smoothing parameter reaches 9.2e+08, the criterion is
genuinely flat, and calling that fit unconverged would be wrong. A
coordinate is reported as sitting at a boundary when its free value
exceeds `edge` AND
[`coord_decrement()`](https://statmodels7.github.io/statmodels7/reference/coord_decrement.md),
what that coordinate alone would still buy, has already met `tol`.
Because of that second condition the threshold cannot change the
verdict, and since 0.127.0 by an exact inequality rather than by an
argument about a maximum: restricting a decrement to a subset of the
coordinates is that same maximum under a constraint and so no larger,
hence a coordinate moved out of the interior set cannot raise what
decides the state. Both `"converged"` and `"boundary"` are certified.
What `edge` decides is how the point is described. The default separates
the measured cases with room on both sides: coordinates that ran to an
edge sit at 9.3, 10.5 and 20.6 on the free scale against 0.13, 0.30 and
2.01 for the ones that did not.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md),
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md),
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)

## Examples

``` r
dd <- data.frame(x = runif(120))
dd$y <- sin(4 * dd$x) + rnorm(120, 0, 0.3)
statmod_certificate(statmod(y ~ s(x, bspline_smooth(k = 8)),
                            distributions7::gaussian1_distrib(), dd))$state
#> [1] "converged"
```
