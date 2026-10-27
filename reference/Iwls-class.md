# The Iterated Weighted Least Squares Method

`Iwls` is the S7 class of the scoring step that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
uses for the jointly fitted smooth block. An object of this class
carries the solver, the budget, the stopping rule and the constants of
the loop. Build one with
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
which validates every argument; this raw constructor does not.

The matrix that the step inverts is not an argument of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md):
it is the argument `inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
stores the matrix that it settled in the properties `hessian`, `backup`,
`switching` and `repair_floor`, which
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
initializes to an unsettled state.

## Usage

``` r
Iwls(
  solver = IwlsSolver(),
  maxit = integer(0),
  tol = integer(0),
  criterion = NULL,
  step_halving = integer(0),
  decrement_limit = integer(0),
  switch_after = integer(0),
  switch_rate = integer(0),
  hessian = distributions7::Information(),
  backup = NULL,
  switching = logical(0),
  repair_floor = integer(0)
)
```

## Arguments

- solver:

  How the step is solved, an object that inherits from
  [IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md):
  [`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
  [`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
  or
  [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md).

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

  The value under which the built-in rule reads the Newton decrement as
  small, a single positive number.

- switch_after:

  A single number: the iteration after which a run on the expected
  information continues on the observed one when the choice was settled
  from
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md).
  `Inf` never switches.

- switch_rate:

  A single number, or `NA`: the ratio of consecutive scores above which
  a run on the expected information counts as contracting slowly. `NA`
  never counts a run as slow.

- hessian:

  The information that the step inverts, once settled: an object that
  inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html),
  set by
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  from its argument `inner_information`. An
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  object stands for an unsettled choice, which
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
  settles against the family.

- backup:

  `NULL`, or the
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  object that takes the step where the observed step finds no acceptable
  point. It is set only by settling an
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  on the observed information.

- switching:

  Whether the run continues on the observed information after
  `switch_after` iterations or a slow contraction. It is set only by
  settling an
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  on the expected information.

- repair_floor:

  The `repair_floor` of the
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  that was settled, or `NA` when the choice was named.

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
iwls()@solver
#> QR solver (tol 1e-07, floor 1e-08)
```
