# The Certificate of the Reported Point

Reads the outer criterion at the reported point, independently of the
path of the search, and returns a state (`"converged"`, `"boundary"`,
`"not converged"` or `"unknown"`) that is determined by the Newton
decrement of the criterion. The result also holds the mode error of the
coefficients, the hyperparameters and flat directions at a boundary, and
the dimensionless readings that
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints.

## Usage

``` r
statmod_certificate(fit, tol = 0.01, flat = 0.002, edge = 8)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- tol:

  The largest rise, in the criterion's own units, that a certified point
  may still have available; see
  [`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md).

- flat:

  The largest curvature \\\lvert A\_{jj}\rvert\\ of the outer criterion
  in a hyperparameter's free value at which that hyperparameter is
  reported as sitting at a boundary; for a coefficient the criterion
  estimates, the curvature times \\\max(1, \gamma_j^2)\\, and for a free
  value of a matrix prior that scales with the matrix, the curvature in
  units of the matrix's scale
  ([`hyper_unit_scale()`](https://statmodels7.github.io/statmodels7/reference/hyper_unit_scale.md)).
  It is also the largest eigenvalue, in absolute value, of a flat
  direction of the equilibrated curvature. It determines only the label
  and never the verdict; see the details.

- edge:

  The free value beyond which a hyperparameter is reported as sitting at
  a boundary where no curvature can be read, and only there.

## Value

A list with `state` (`"converged"`, `"boundary"`, `"not converged"` or
`"unknown"`), `decrement`, `gradient`, `mode_error`, `curvature`,
`boundary`, `boundary_key`, `ridge` and `reason`. `ridge` is a list with
one named vector per flat direction of the outer curvature, over the
coordinates that move along it. The verdict uses `decrement`, and
`gradient` is reported beside it. `curvature` is `"analytic"` when the
decrement was read against an analytic Hessian, `"differenced"` against
one differenced from the exact gradient, `"differenced criterion"` when
the criterion itself was differenced by refits, and `NA` when no
curvature was read. `boundary_key` names the hyperparameters among the
coordinates of `boundary`, in the key by which
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
labels its rows;
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reads it to leave a standard error off a hyperparameter pinned at an
edge. `boundary_coef` names the coefficients estimated on the outer
criterion among the coordinates of `boundary`, as `parameter:name`, and
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reads it to leave their standard errors, tests and intervals off. The
state is `"not converged"` also when the fit recorded a coordinate at
the edge of its chart whose return raises the log-likelihood.

Three further fields are those that
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints, and its page describes how they are read: `inner` and `outer`,
each a named vector with `gradient` and `eigen` from
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
(`NA` where there is nothing to read), and `zeros`, which reads `NULL`
where no coefficient sits at the kink of a penalty and is otherwise a
named vector with `ratio` and `n` from
[`zero_readings()`](https://statmodels7.github.io/statmodels7/reference/zero_readings.md).

## Details

The flag of an optimizer records whether a search stopped on its own
rule, which describes the search and does not order fits by quality. The
certificate describes the reported point and does not depend on the path
that led to it.

The state comes from the rise that the criterion would still buy, and
the mode error is reported beside it. That rise is the Newton decrement
\\\tfrac{1}{2}g^\top A^{-1}g\\ of
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md),
in the criterion's own units, and `tol` is a threshold on it. The mode
error is the number of log-likelihood units by which the coefficients
lie above the penalized mode. Where an outer criterion exists, the mode
error does not enter the verdict; for a model with no penalty and no
coefficient estimated on a criterion there is no other reading, and the
mode error against
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
determines the state (see below).

The verdict is a threshold on the Newton decrement and not on the
gradient, because the gradient of a criterion summed over \\n\\
observations and \\p\\ penalized coefficients depends on the model, so
that a single threshold cannot serve every model. See
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md).
The default `tol` is 1e-2. A user can see and check a report of
`"not converged"` at a good point, whereas a user cannot check the
certification of a point that is not a maximum, so the threshold is set
low.

The cost of the certificate is one outer gradient, one curvature and one
solve, at a point that the fit already holds. The criterion is
reconstructed from `fit@spec`, so the reading refers to the fitted
model. Where the form carries an analytic outer Hessian, the certificate
does not refit. Where it does not (a criterion on the expected
information for a family that does not supply the second derivative of
its expected information, a model with a refreshable block or a
structural term, or a penalty whose Hessian moves with the coefficients
beside a structural term), the curvature is obtained by central
differences of the exact gradient over the hyperparameters and the
coefficients that the criterion estimates, which takes \\4(n_h +
n\_\gamma)\\ refits.
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md)
records which route was taken, and the result reports it in `curvature`.

Where there is no outer gradient to read, the state depends on the case.
A model with no penalty and no coefficient that the criterion estimates
(a family of one parameter, or `marginal` empty), as with `linpar`,
`nl`, `seg`, `jump` and `jseg`, has no coordinate for which an outer
gradient could be read, so the only question left is whether the inner
fit reached its mode. The state is `"converged"` when the mode error is
at most
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
and `"not converged"` otherwise. In every other case without an outer
gradient the state is `"unknown"` and `reason` names the case: a model
whose only hyperparameters are kinked (`lasso`, `scad`, `mcp`), which
are chosen along a path because a Laplace approximation at a mode that
sits on the kink has no meaning, so that the hyperparameter is the
argument of a minimum over a grid; hyperparameters that are all held;
hyperparameters chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md); and
a model fitted with no outer criterion. `"unknown"` states that the
certificate has no outer gradient to read and is not a report of failure
to converge. With a kinked penalty the mode error is reported over the
coordinates that the kink leaves free
([`free_of_kinks()`](https://statmodels7.github.io/statmodels7/reference/free_of_kinks.md)):
at a coefficient that the penalty has set to zero the score does not
vanish but lies in the subdifferential, so that reading the score there
would give a false mode error. A model that also carries a smooth
penalty estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
has an outer gradient, which the certificate reads as usual.

A form whose criterion has no exact gradient
([`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
at order 1) has the criterion itself differenced by refits when it has
at most two hyperparameters and no coefficient estimated on the
criterion (three refits for one coordinate, nine for two), and
`curvature` is then `"differenced criterion"`. Otherwise the state is
`"unknown"`, and `reason` gives the cause. This is consistent with
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md),
which differences a Hessian where the analytic one is missing, because
that function differences the exact gradient.

A hyperparameter may run to an edge and belong there: on a covariate
with no effect the smoothing parameter grows without bound, the
criterion is flat in it, and the fit has not failed to converge. A
coordinate is reported as sitting at a boundary where the criterion no
longer moves with it, that is where its own curvature \\\lvert
A\_{jj}\rvert\\ is at most `flat`. Such a hyperparameter has no
meaningful interval, because \\1/\sqrt{A\_{jj}}\\ is the conditional
standard error of its free value, and `boundary_key` gives
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
the keys of these hyperparameters, so that it leaves their standard
errors out.

The label is read from the curvature and not from the free value,
because the size of a smoothing parameter depends on how its penalty is
normalized and on the units of the response. The default `flat` of 2e-3
corresponds to a conditional standard error of the free value above 22.

The set of coordinates excluded from the verdict is narrower: a
coordinate belongs to it exactly when it is reported at a boundary and
its
[`coord_decrement()`](https://statmodels7.github.io/statmodels7/reference/coord_decrement.md),
the rise that the coordinate would still give if it were the only one to
move, is at most `tol`. The second condition prevents `flat` from
determining the state: restricting a decrement to a subset of the
coordinates is that same maximum under a constraint and so no larger,
hence a coordinate moved out of the interior set cannot raise the
quantity that determines the state. A coordinate at an edge that has not
settled stays under test and is also named. Both `"converged"` and
`"boundary"` are certified.

The set of reported coordinates and the set excluded from the verdict
can differ: at a boundary the curvature collapses along with the
gradient, so the one-coordinate reading \\g_j^2/(2A\_{jj})\\ is a ratio
of two quantities going to zero and can be large. Such a coordinate is
named and stays under test.

A coefficient that the criterion estimates beside the hyperparameters
(the argument `marginal` of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md))
is read the same way, but its curvature is multiplied by \\\max(1,
\gamma_j^2)\\: it is reported at a boundary where its conditional
standard error exceeds \\22\\ times \\\max(1, \lvert\gamma_j\rvert)\\. A
coefficient carries the units of its parameter and of its covariate,
which a hyperparameter's free value does not, so the bare curvature
depends on those units. The ratio of the coefficient to its conditional
standard error does not. The degrees of freedom of a Student \\t\\ at
their gaussian limit are an example of such a coordinate.

Where no curvature can be read at all (both routes fail to produce one,
and the state is then `"unknown"`) the label falls back on the free
value, `edge`.

A ridge is a direction along which several coordinates move together
while the criterion does not change, although the criterion is curved in
each of those coordinates separately. It is read on the equilibrated
curvature \\A/(ss')\\, \\s_j = \sqrt{\|A\_{jj}\|}\\, whose eigenvalues
are dimensionless. An eigenvalue of at most `flat` in absolute value
names its eigenvector, written on the original coordinates in `boundary`
and in `ridge` (see
[`certificate_ridges()`](https://statmodels7.github.io/statmodels7/reference/certificate_ridges.md)).
A ridge whose own decrement is within `tol` leaves the verdict and the
state is `"boundary"`. A Burr XII shape that runs to its Weibull limit
in one level of a factor is an example: the eigenvector moves only that
level.

A fit that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
has recorded with a coordinate at the edge of its chart whose return
raises the log-likelihood has the state `"not converged"` whatever the
readings above give, and `reason` names the coordinate.

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
