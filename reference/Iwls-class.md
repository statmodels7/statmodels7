# The Iterated Weighted Least Squares Method

`Iwls` is the S7 class of the scoring step that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
uses for the jointly fitted smooth block. An object of this class
carries the whole configuration of that step: which curvature to invert,
how to solve the system, how long to run and when to stop. Build one
with
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
which validates every argument; this raw constructor does not.

## Usage

``` r
Iwls(
  hessian = character(0),
  approx = character(0),
  decomposition = character(0),
  maxit = integer(0),
  tol = integer(0),
  criterion = NULL,
  step_halving = integer(0),
  fallback = logical(0),
  switch_after = integer(0)
)
```

## Arguments

- hessian:

  `"expected"` for Fisher scoring, `"observed"` for Newton, or `"auto"`,
  which
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  turns into one of the two once the distribution is known. The value is
  a single string.

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
  tries only the full step.

- fallback:

  A single logical: whether the expected information takes the step
  where the observed step finds no acceptable point. An observed
  penalized information that is not positive definite is first repaired
  by flooring its eigenvalues, and the step on the repaired matrix is
  tried.
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  sets this argument to `TRUE` for `"auto"` on a family that it does not
  settle on the expected information.

- switch_after:

  A single number: the iteration after which a run on the expected
  information continues on the observed one, unless it switched earlier
  because the score contracted slowly (see
  [`iwls_slow()`](https://statmodels7.github.io/statmodels7/reference/iwls_slow.md)).
  The expected information then takes the step only where the observed
  step finds no acceptable point. `Inf` never switches.
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  sets it, for `"auto"` on a family whose expected information is exact
  and not costly to evaluate, to
  [`iwls_switch_after()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_after.md).

## Value

An object of class `Iwls` with one property per argument above, and each
property holds the value that was passed.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
the constructor to use.

## Examples

``` r
S7::S7_inherits(iwls(), Iwls)
#> [1] TRUE
iwls()@decomposition
#> [1] "qr"
```
