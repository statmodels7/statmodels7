# Iterated Weighted Least Squares

The scoring step, written out, with the curvature and the decomposition
left to the caller.

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

  Either `"expected"` – Fisher scoring – or `"observed"`, which is
  Newton.

- approx:

  The approximation of the expected information, read only where the
  family has no closed form.

- decomposition:

  How the step is solved: `"qr"`, `"svd"`, `"chol"` or
  `"chol_crossprod"`.

- maxit:

  The iteration budget.

- tol:

  The stopping tolerance on the scaled score.

- criterion:

  An optimizers7 `criterion` driving the loop in place of `tol`, or
  `NULL` for the built-in rule.

- step_halving:

  The number of halvings allowed before a step is abandoned.

- x:

  An
  [`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- ...:

  Unused.

## Value

An object of class
[`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md).

## Details

The name is `iwls` and not `irls`: the re-weighting is already implied
by the iteration, so the second `r` says nothing that *iterated* has not
said.

**The curvature.** `hessian = "expected"` inverts the expected
information, which is Fisher scoring; `"observed"` inverts the observed
Hessian, which is Newton. `approx` is handed to distributions7 and is
read only where the family has no closed expected information; asking
for one where it would be ignored is refused, since an argument accepted
and ignored is worse than one that errors.

**The decomposition.** The step solves \\(X'WX + S)\delta = X'g -
S\beta\\, and how depends on this argument.

- `"qr"`:

  a QR of the augmented design \\\[\sqrt{W}X;\\ \mathrm{chol}(S)\]\\.
  The default, and not as a matter of taste: forming \\X'X\\ squares a
  conditioning that the break-point terms bound to \\\epsilon^{-1/2}\\
  by construction, and it is why a closed-form ridge is beaten at \\p =
  200\\ by a method that never forms it.

- `"svd"`:

  a singular value decomposition of the same augmented design, which
  reports the numerical rank rather than failing on a deficient block.

- `"chol"`:

  a Cholesky factor of the penalized information, assembled without the
  augmentation.

- `"chol_crossprod"`:

  the same, forming \\X'WX\\ explicitly. The fastest per iteration and
  the worst conditioned; offered because the choice is the user's.

**The stopping rule.** `iwls` is a scoring step and not an optimizer, so
it carries its own loop, but the rule that ends it is the caller's to
choose. With `criterion = NULL` the built-in rule applies: the score per
observation \\\max_j\lvert g_j\rvert / n\\ against `tol`, with the
dimensionless reading of
[`iwls_score`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
arbitrating the final verdict. Any optimizers7 `criterion` may drive the
loop instead, and then `tol` is not read at all, so passing both is an
error rather than a silent choice between them.

What the rule is shown is the state
[`optimizers7::crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.html)
documents, with two things worth knowing. Its `gradient` is the score
PER OBSERVATION, the quantity the built-in rule compares, so
`criterion = optimizers7::crit_grad(t)` is `tol = t` exactly and a
threshold means the same at \\n = 10\\ and at \\n = 10^7\\. Its
objective is the penalized log-likelihood UNAVERAGED, which is the scale
the penalty is added on, so a rule reading the objective's absolute
value rather than its relative change carries the sample size with it.
On the first iteration `f_old` and `x_old` are `NULL`, there being no
previous point, so a rule reading a change returns `FALSE` there by
construction. A rule needing something the step does not compute – a
stationarity measure, which belongs to the derivative-free methods – is
rejected at construction rather than sitting there never firing.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls_score`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)

## Examples

``` r
iwls()
#> iwls: expected information, qr
#>   maxit 100, tol 1e-06
iwls(hessian = "observed", decomposition = "svd")
#> iwls: observed information, svd
#>   maxit 100, tol 1e-06
iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
                                       optimizers7::crit_rel_obj(1e-12)))
#> iwls: expected information, qr
#>   maxit 100, gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
```
