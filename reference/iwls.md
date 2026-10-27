# Iterated Weighted Least Squares

Configures the scoring step with which
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits the smooth block, and returns it as an object to pass as
`inner_optimizer`. This is the default, so `statmod(...)` and
`statmod(..., inner_optimizer = iwls())` fit the same model; call it by
name to change the solver, the budget, the stopping rule or the
constants of the loop.

The matrix that the step inverts is chosen by the argument
`inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
(see
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)),
and not by `iwls()`: the step is a Newton method on whatever curvature
the objective offers. Every argument of `iwls()` has a default, and the
object carries no data, so one `iwls()` can be passed to any number of
fits.

## Usage

``` r
iwls(
  solver = solve_qr(),
  maxit = 1000L,
  tol = 1e-06,
  criterion = NULL,
  step_halving = 30L,
  decrement_limit = 1e-06,
  switch_after = 10,
  switch_rate = 0.5
)

# S3 method for class 'Iwls'
print(x, ...)
```

## Arguments

- solver:

  The factorization with which each step is solved, an object that
  inherits from
  [IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md):
  [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
  (the default),
  [`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
  or
  [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md).
  [iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
  compares the three and gives their tolerances.

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
  tries only the full step.

- decrement_limit:

  The value under which the Newton decrement is read as small by the
  built-in stopping rule, a single positive number. The default is
  `1e-6`.

- switch_after:

  The iteration after which a run settled on the expected information
  from
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  continues on the observed information, a single positive number; `Inf`
  never switches. The default is `10`.

- switch_rate:

  The ratio of consecutive scores above which two consecutive steps
  count as a slow contraction of Fisher scoring, a single number or
  `NA`, which never counts a run as slow. The default is `0.5`.

- x:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- ...:

  Unused.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object holding the settings given and no data. Its `hessian` is an
unsettled
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
until
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
settles the argument `inner_information`.

## The name

The name abbreviates iterated weighted least squares. The word iterated
already implies the re-weighting, so the abbreviation `irls` is not
used.

## The curvature

The step inverts the matrix that `inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
selects. With
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
it inverts the expected information, which is Fisher scoring, and with
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
the observed Hessian, which is Newton. The expected information is
positive definite wherever the family is regular, so Fisher scoring
takes usable steps from a poor start where Newton's matrix may be
indefinite; Newton converges faster near the optimum. An explicit choice
is the only matrix that the step uses: the observed information is never
replaced by the expected one.

The default `inner_information = auto_info()` is settled against the
family when the fit starts. The expected information is used where the
family writes it exactly and at an acceptable cost. Otherwise the
observed information is used, and the expected information takes the
step where the observed step, after the repair of an indefinite matrix,
finds no acceptable point. A trial point whose objective raises an error
is read as rejected. The fit records the settled choice, so
`fit@methods$smooth` shows which curvature ran.

On a family that
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
settles on the expected information, a run continues on the observed
information as soon as Fisher scoring contracts the score slowly, which
means that each of two consecutive ratios of the score is above
`switch_rate`, and at the latest after `switch_after` scoring steps.
Fisher scoring converges only linearly near the mode wherever the two
informations differ there, which they do whenever the dispersion has an
equation of its own, whereas Newton's step converges quadratically.

Where the stopping rule is met, the run takes one more full step. The
step uses the observed information if the choice was settled on the
expected information from
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
and has not yet switched, and otherwise the information that the run is
using. The marginal criterion carries \\\log\lvert K(\beta)\rvert\\,
which is not stationary in \\\beta\\, so a mode left short by \\\delta\\
enters it at first order; a warm-started fit whose score is already
under the rule would stop without moving and leave the criterion out by
an amount linear in the outer step. The extra step takes the error to
second order. It is kept unless it raises the objective by more than a
relative \\10^{-8}\\: a sufficient-decrease test at a point already
under the rule would accept or reject by the last bits, and the
criterion would jump between nearby hyperparameters.

The `method` carried by
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
and
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
is passed to distributions7. With `method = NULL`, the default, the
closed form of the family is used where the family has one
([`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
records which case a family is in), and
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
the outer product of the observed scores, which costs one gradient,
otherwise. An object is used for every family, also for one that has the
closed form.
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
and
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
evaluate the expectation itself (a sum over the support or a quadrature,
or simulated draws) and cost much more inside a loop that rebuilds the
curvature at every iteration.
[distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)
compares the three. The price of the cheap matrix is paid in variance
and not in the answer: the score is exact, so a scoring step driven by
any positive definite matrix reaches the same stationary point, and the
outer product is positive semidefinite by construction, whereas the
observed Hessian need not be. The information behind a standard error is
chosen separately; see
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

## The solver

The step solves \\(X'WX + S)\delta = X'g - S\beta\\, and the argument
`solver` selects the method.

- [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md):

  a QR of the augmented design \\\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\,
  and the default. Forming \\X'X\\ squares the condition number, and the
  design of a break-point term is bounded only to \\\epsilon^{-1/2}\\ by
  construction, so that squaring it exhausts double precision.

- [`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md):

  a singular value decomposition of the same augmented design. It
  reports a numerical rank where a deficient block would make a Cholesky
  fail.

- [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md):

  a Cholesky factor of the penalized information \\X'WX + S\\, assembled
  without the augmentation. The information is formed explicitly, so the
  condition number is squared.

[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
and
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
need a square root of the per-observation curvature and of the penalty.
Where either does not exist, the step uses the assembled information, as
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
does.
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
compares the three and gives the tolerances that each one takes.

## The stopping rule

A scoring step is not an optimizer and carries its own loop, but the
rule that ends the loop belongs to the caller.

With `criterion = NULL`, the default, the built-in rule applies. The
rule is met when the score per observation \\\max_j \lvert g_j \rvert /
n\\, taken over the free coordinates, is below `tol` and the Newton
decrement \\g^\top (H+S)^{-1} g / 2\\ is at most `decrement_limit`. The
decrement is in log-likelihood units and does not depend on the scale of
the response. When the loop ends without meeting the rule, the
dimensionless reading of
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
can relabel the stopped run as converged. Any optimizers7 criterion may
drive the loop instead. It replaces the whole built-in rule, including
the decrement condition and the dimensionless reading, and `tol` is then
not read, so passing both is an error.

The state that a criterion receives is the one that
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

- On the first iteration `f_old` and `x_old` are `NULL`, because there
  is no previous point, so a rule reading a change returns `FALSE` there
  by construction.

A criterion that needs a quantity that the step does not compute is
rejected at construction. A stationarity measure is an example: it
belongs to the derivative-free methods, and a rule that requires one
would remain in the loop without ever being met.

## The constants of the loop

Three constants of the loop are arguments. `decrement_limit` is the
value under which the built-in rule reads the Newton decrement as small
(default `1e-6`, three orders under the reading that the outer search
uses to decide whether a criterion is usable). `switch_after` is the
iteration after which a run settled on the expected information
continues on the observed one (default `10`; `Inf` never switches).
`switch_rate` is the ratio of consecutive scores above which two
consecutive steps count as a slow contraction (default `0.5`; `NA` never
counts one). The last two are read only when `inner_information` is
settled from
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
on the expected information.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the argument `inner_information` and for where the object is passed,
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
and
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the matrix and the solver,
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
for the built-in rule's final verdict,
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
for the criteria that can replace it.

## Examples

``` r
iwls()
#> iwls: qr
#>   maxit 1000, tol 1e-06

# This call solves through an SVD, which reports a numerical rank where a
# deficient block would make a Cholesky fail.
iwls(solver = solve_svd())
#> iwls: svd
#>   maxit 1000, tol 1e-06

# This call uses an optimizers7 stopping rule in place of tol. The
# gradient in the state that the rule receives is the score per
# observation, as for tol.
iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
                                       optimizers7::crit_rel_obj(1e-12)))
#> iwls: qr
#>   maxit 1000, gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)

# The matrix is chosen in statmod(): this is Newton, on the observed
# information.
set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_optimizer = iwls(solver = solve_svd()),
               inner_information = distributions7::observed_info())
fit@methods$smooth
#> iwls: svd, observed information
#>   maxit 1000, tol 1e-06

# Both at once is an error: tol would not be read.
try(iwls(tol = 1e-8, criterion = optimizers7::crit_grad(1e-8)))
#> Error : 'tol' and 'criterion' both set the stopping rule: pass one.
#>   'tol' belongs to the built-in rule, which a criterion replaces.
```
