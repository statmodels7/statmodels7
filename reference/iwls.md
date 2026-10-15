# Iterated Weighted Least Squares

Configures the scoring step with which
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits the smooth block, and returns it as an object to pass as
`inner_optimizer`. This is the default, so `statmod(...)` and
`statmod(..., inner_optimizer = iwls())` fit the same model; call it by
name to change the curvature, the decomposition, the budget or the
stopping rule.

Every argument has a default, and the object carries no data: one
`iwls()` can be passed to any number of fits.

## Usage

``` r
iwls(
  hessian = c("auto", "expected", "observed"),
  approx = c("opg", "bartlett", "integrate", "mc"),
  decomposition = c("qr", "svd", "chol"),
  maxit = 1000L,
  tol = 1e-06,
  criterion = NULL,
  step_halving = 30L
)

# S3 method for class 'Iwls'
print(x, ...)
```

## Arguments

- hessian:

  `"expected"` for Fisher scoring, `"observed"` for Newton, or `"auto"`,
  which
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  turns into one of the two once the distribution is known. A single
  string.

- approx:

  How the expected information is approximated for a family with no
  closed form: `"opg"`, `"bartlett"`, `"integrate"` or `"mc"`.

- decomposition:

  How the step is solved: `"qr"`, `"svd"` or `"chol"`.

- maxit:

  The iteration budget, a single positive whole number.

- tol:

  The stopping tolerance on the score per observation, a single positive
  number.

- criterion:

  An optimizers7 `criterion` object driving the loop in place of `tol`,
  or `NULL` for the built-in rule.

- step_halving:

  The number of halvings tried after the full step before a step is
  abandoned, a single non-negative whole number. The line search tries
  the full step and then at most `step_halving` halved steps, so `0`
  tries the full step alone.

- x:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- ...:

  Unused.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object holding the settings given and no data. Its `fallback` is `FALSE`
and its `switch_after` is `Inf` until
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
settles an `"auto"`.

## The name

The name abbreviates iterated weighted least squares. The word iterated
already implies the re-weighting, so the abbreviation `irls` is not
used.

## The curvature

`hessian = "expected"` inverts the expected information, which is Fisher
scoring, and `"observed"` inverts the observed Hessian, which is Newton.
The expected information is positive definite wherever the family is
regular, so Fisher scoring takes usable steps from a poor start where
Newton's matrix may be indefinite; Newton converges faster near the
optimum and is what the exact outer gradient needs.

`"auto"`, the default, is settled against the family when the fit
starts, by
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md).
The expected information is used where the family writes it exactly and
at an acceptable cost. Otherwise the observed information is used. An
observed penalized information that is not positive definite is first
repaired by flooring its eigenvalues, and the expected information takes
the step only where the repaired step finds no acceptable point. A trial
point whose objective raises an error is read as rejected. The fit
records the settled method, so `fit@methods$smooth` shows which
curvature ran.

On a family that
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
settles on the expected information from `"auto"`, a run continues on
the observed information as soon as Fisher scoring contracts the score
slowly, each of two consecutive ratios of the score being above
[`iwls_switch_rate()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_rate.md),
and at the latest after
[`iwls_switch_after()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_after.md)
scoring steps. The expected information then takes the step only where
the observed step finds no acceptable point. Fisher scoring converges
only linearly near the mode wherever the two informations differ there,
which they do whenever the dispersion has an equation of its own,
whereas Newton's step converges quadratically.

Where the rule is met, the run takes one more full step. The step uses
the observed information if the method was settled on the expected
information from `"auto"` and has not yet switched, and otherwise the
information that the run is using. The marginal criterion carries
\\\log\lvert K(\beta)\rvert\\, which is not stationary in \\\beta\\, so
a mode left short by \\\delta\\ enters it at first order; a warm-started
fit whose score is already under the rule would stop without moving and
leave the criterion out by an amount linear in the outer step. The extra
step takes the error to second order. It is kept unless it raises the
objective by more than a relative \\10^{-8}\\: a sufficient-decrease
test at a point already under the rule would accept or reject by the
last bits, and the criterion would jump between nearby hyperparameters.

`approx` reaches distributions7 and is read only where the family has no
closed expected information; elsewhere the family's own method answers
and the argument does not enter. The argument is accepted for every
family because an `iwls()` object is built before it meets a
distribution and can be passed to any number of fits, so whether the
family has a closed expected information cannot be checked at
construction.
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is the predicate that records which case a family is in.

Its default is `"opg"`, the outer product of the observed scores, which
costs one gradient. The alternative readings evaluate the expectation
itself (a sum over the support for a discrete family, a quadrature for a
continuous one) and cost much more inside a loop that rebuilds the
curvature at every iteration. What the cheap matrix gives up is variance
and not the answer: the score is exact, so a scoring step driven by any
positive definite matrix reaches the same stationary point, and the
outer product is positive semidefinite by construction where the
observed Hessian need not be. The information behind a standard error is
chosen separately; see
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## The decomposition

The step solves \\(X'WX + S)\delta = X'g - S\beta\\, and this argument
selects the method.

- `"qr"`:

  a QR of the augmented design \\\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\,
  and the default. Forming \\X'X\\ squares the condition number, and the
  design of a break-point term is bounded only to \\\epsilon^{-1/2}\\ by
  construction, so that squaring it exhausts double precision.

- `"svd"`:

  a singular value decomposition of the same augmented design. It
  reports a numerical rank where a deficient block would make a Cholesky
  fail.

- `"chol"`:

  a Cholesky factor of the penalized information \\X'WX + S\\, assembled
  without the augmentation. The information is formed explicitly, so the
  condition number is squared.

The `"qr"` and `"svd"` methods need a square root of the per-observation
curvature and of the penalty. Where either does not exist, the step uses
the assembled information, as `"chol"` does.

## The stopping rule

A scoring step is not an optimizer and carries its own loop, but the
rule that ends the loop belongs to the caller.

With `criterion = NULL`, the default, the built-in rule applies. The
rule is met when the score per observation \\\max_j \lvert g_j \rvert /
n\\, taken over the free coordinates, is below `tol` and the Newton
decrement \\g^\top (H+S)^{-1} g / 2\\ is at most \\10^{-6}\\. The
decrement is in log-likelihood units and does not depend on the scale of
the response. When the loop ends without meeting the rule, the
dimensionless reading of
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
can relabel the stopped run as converged. Any optimizers7 criterion may
drive the loop instead. It replaces the whole built-in rule, the
decrement condition and the dimensionless reading included, and `tol` is
then not read, so passing both is an error.

The state that a criterion receives is the one
[`optimizers7::crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.html)
documents. In this loop the state has three properties:

- Its `gradient` is the score per observation, the quantity that the
  built-in rule compares with `tol`. The threshold of
  `optimizers7::crit_grad(t)` therefore applies to the same quantity as
  `tol = t` and has the same meaning at \\n = 10\\ and at \\n = 10^7\\.

- Its `f_new` is the penalized objective, the negative penalized
  log-likelihood, not divided by the sample size, which is the scale on
  which the penalty is added. A rule that reads the absolute value of
  the objective therefore depends on the sample size; a rule that reads
  its relative change does not.

- On the first iteration `f_old` and `x_old` are `NULL`, there being no
  previous point, so a rule reading a change returns `FALSE` there by
  construction.

A criterion that needs a quantity the step does not compute is rejected
at construction. A stationarity measure is an example: it belongs to the
derivative-free methods, and a rule that requires one would remain in
the loop without ever being met.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where it is passed,
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
for the built-in rule's final verdict,
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
for the criteria that can replace it.

## Examples

``` r
iwls()
#> iwls: auto information, qr
#>   maxit 1000, tol 1e-06

# Newton instead of Fisher scoring, solved through an SVD, which reports a
# numerical rank where a deficient block would make a Cholesky fail.
iwls(hessian = "observed", decomposition = "svd")
#> iwls: observed information, svd
#>   maxit 1000, tol 1e-06

# An optimizers7 stopping rule in place of tol. The gradient in the state
# that the rule receives is the score per observation, as for tol.
iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
                                       optimizers7::crit_rel_obj(1e-12)))
#> iwls: auto information, qr
#>   maxit 1000, gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)

# Both at once is an error: tol would not be read.
try(iwls(tol = 1e-8, criterion = optimizers7::crit_grad(1e-8)))
#> Error : 'tol' and 'criterion' both say when to stop: pass one.
#>   'tol' belongs to the built-in rule, which a criterion replaces.
```
