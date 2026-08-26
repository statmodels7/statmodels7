# The Iterated Weighted Least Squares Method

The S7 class of the scoring step
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
  step_halving = integer(0)
)
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

## Value

An object of class `Iwls` with one property per argument above, each
holding what was passed.

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
