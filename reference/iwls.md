# Iterated Weighted Least Squares

Configures the scoring step
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits the smooth block with, and returns it as an object to pass as
`inner_optimizer`. This is the default, so `statmod(...)` and
`statmod(..., inner_optimizer = iwls())` fit the same model; call it by
name to change the curvature, the decomposition, the budget or the
stopping rule.

Every argument has a default, and the object carries no data: one
`iwls()` can be passed to any number of fits.

## Usage

``` r
iwls(
  hessian = c("expected", "observed"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  decomposition = c("qr", "svd", "chol", "chol_crossprod"),
  maxit = 100L,
  tol = 1e-06,
  criterion = NULL,
  step_halving = 30L
)

# S3 method for class 'Iwls'
print(x, ...)
```

## Arguments

- hessian:

  `"expected"` for Fisher scoring or `"observed"` for Newton, a single
  string.

- approx:

  How the expected information is approximated for a family with no
  closed form: `"bartlett"`, `"integrate"` or `"mc"`.

- decomposition:

  How the step is solved: `"qr"`, `"svd"`, `"chol"` or
  `"chol_crossprod"`.

- maxit:

  The iteration budget, a single positive number.

- tol:

  The stopping tolerance on the score per observation.

- criterion:

  An optimizers7 `criterion` object driving the loop in place of `tol`,
  or `NULL` for the built-in rule.

- step_halving:

  How many halvings are allowed before a step is abandoned, a single
  non-negative number.

- x:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- ...:

  Unused.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object, holding the seven settings above and no data.

## The name

`iwls`, not `irls`. The iteration already implies the re-weighting, so
the second `r` says nothing *iterated* has not said.

## The curvature

`hessian = "expected"` inverts the expected information, which is Fisher
scoring, and `"observed"` inverts the observed Hessian, which is Newton.
The expected information is positive definite wherever the family is
regular, so Fisher scoring takes usable steps from a poor start where
Newton's matrix may be indefinite; Newton converges faster near the
optimum and is what the exact outer gradient needs.

`approx` reaches distributions7 and is read only where the family has no
closed expected information. Asking for one where it would be ignored is
an error: an argument accepted and ignored is worse than one that
refuses.

## The decomposition

The step solves \\(X'WX + S)\delta = X'g - S\beta\\, and this argument
says how.

- `"qr"`:

  a QR of the augmented design \\\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\,
  and the default. Forming \\X'X\\ squares the condition number, and a
  break-point term's design is bounded only to \\\epsilon^{-1/2}\\ by
  construction, so squaring it exhausts double precision. This is also
  why a closed-form ridge is beaten at \\p = 200\\ by a method that
  never forms \\X'X\\.

- `"svd"`:

  a singular value decomposition of the same augmented design. It
  reports a numerical rank where a deficient block would make a Cholesky
  fail.

- `"chol"`:

  a Cholesky factor of the penalized information, assembled without the
  augmentation.

- `"chol_crossprod"`:

  the same, forming \\X'WX\\ explicitly. The fastest per iteration and
  the worst conditioned; offered because the choice is the user's.

## The stopping rule

A scoring step is not an optimizer and carries its own loop, but the
rule that ends the loop belongs to the caller.

With `criterion = NULL`, the default, the built-in rule applies: the
score per observation \\\max_j \lvert g_j \rvert / n\\ against `tol`,
with the dimensionless reading of
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
arbitrating the final verdict. Any optimizers7 criterion may drive the
loop instead, and then `tol` is never read, so passing both is an error.

What a criterion is shown is the state
[`optimizers7::crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.html)
documents, with three things worth knowing about it here:

- Its `gradient` is the score **per observation**, the quantity the
  built-in rule compares. So `criterion = optimizers7::crit_grad(t)` is
  exactly `tol = t`, and either threshold means the same at \\n = 10\\
  and at \\n = 10^7\\.

- Its objective is the penalized log-likelihood **unaveraged**, which is
  the scale the penalty is added on. A rule reading the objective's
  absolute value therefore carries the sample size with it; one reading
  its relative change does not.

- On the first iteration `f_old` and `x_old` are `NULL`, there being no
  previous point, so a rule reading a change returns `FALSE` there by
  construction.

A criterion needing something the step does not compute is rejected at
construction. A stationarity measure is the case: it belongs to the
derivative-free methods, and a rule asking for one would otherwise sit
in the loop and never fire.

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
#> iwls: expected information, qr
#>   maxit 100, tol 1e-06

# Newton instead of Fisher scoring, solved through an SVD, which reports a
# numerical rank where a deficient block would make a Cholesky fail.
iwls(hessian = "observed", decomposition = "svd")
#> iwls: observed information, svd
#>   maxit 100, tol 1e-06

# An optimizers7 stopping rule in place of tol. crit_grad(t) here is
# exactly tol = t, the state's gradient being the score per observation.
iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
                                       optimizers7::crit_rel_obj(1e-12)))
#> iwls: expected information, qr
#>   maxit 100, gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)

# Both at once is an error: tol would be read by nobody.
try(iwls(tol = 1e-8, criterion = optimizers7::crit_grad(1e-8)))
#> Error : 'tol' and 'criterion' both say when to stop: pass one.
#>   The rule reads the score per observation, so crit_grad(tol) is what 'tol' means.
```
