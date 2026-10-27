# A Test of One Coefficient Against One Value

Holds the result of testing that one coefficient of a fit equals one
value: the statistic, its distribution, the p-value and what was tested.

## Usage

``` r
StatmodTest(
  test = test_statistic(),
  statistic = integer(0),
  df = integer(0),
  p.value = integer(0),
  estimate = integer(0),
  null_value = integer(0),
  param = character(0),
  coefficient = character(0),
  mode_error = integer(0),
  converged = logical(0)
)
```

## Arguments

- test:

  The
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
  object that was computed:
  [`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
  [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
  [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
  or
  [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md).
  Its `label` is the name under which
  [`print()`](https://rdrr.io/r/base/print.html) shows the statistic.

- statistic:

  Its value, compared with a \\\chi^2_1\\.

- df:

  One, always.

- p.value:

  The upper tail of that \\\chi^2_1\\.

- estimate:

  The unrestricted estimate of the coefficient.

- null_value:

  The value against which the coefficient was tested.

- param, coefficient:

  Which coefficient of which equation.

- mode_error:

  How far above its own mode the restricted refit stopped, in
  log-likelihood units; `NA` for the Wald statistic.

- converged:

  The inner optimizer's flag at the restricted refit; `NA` for the Wald
  statistic.

## Value

An object of class `StatmodTest`.

## Details

The four statistics are compared with the same \\\chi^2_1\\, so `df` is
always 1 and the object holds one kind of quantity whatever produced it.

## Whether the refit reached its mode

`mode_error` and `converged` both describe the restricted refit and both
are `NA` for the Wald statistic, which reads the unrestricted fit and
refits nothing. They measure different things, and only the first is
printed.

`converged` is the inner optimizer's flag: whether a stopping rule
fired, a boolean about a threshold on a score whose size depends on the
model. `mode_error` is
[`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
the decrease of the penalized objective that a Newton step from the
point returned by the refit would predict, in log-likelihood units,
judged against
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).
The two need not agree, because a stopping rule reads a change in the
objective whereas the mode error reads the distance to the mode. See
[`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md).

A statistic computed from a refit that stopped short is the value at the
point where the refit stopped.
[`print()`](https://rdrr.io/r/base/print.html) reports this, reading
`mode_error`.

The class declares its properties, as every object that this toolkit
returns does, whereas the `htest` returned by the tests of base R is a
list with a class attribute and no declared structure, so that it cannot
be validated.

## See also

[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
which builds one, and
[`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
which `mode_error` reports.

## Examples

``` r
dd <- data.frame(x = seq(-1, 1, length.out = 40))
dd$y <- rpois(40, exp(0.3 + dd$x))
fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
S7::S7_inherits(statmod_test(fit, param = "mu", coefname = "x", value = 1),
                StatmodTest)
#> [1] TRUE
```
