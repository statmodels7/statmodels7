# A Test of One Coefficient Against One Value

The result of testing that one coefficient of a fit equals one value,
carrying the statistic, its distribution, the p-value and what was
tested.

## Usage

``` r
StatmodTest(
  test = character(0),
  statistic = integer(0),
  df = integer(0),
  p.value = integer(0),
  estimate = integer(0),
  null_value = integer(0),
  parameter = character(0),
  coefficient = character(0),
  mode_error = integer(0),
  converged = logical(0)
)
```

## Arguments

- test:

  Which of the four statistics, as
  [`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
  was asked for it.

- statistic:

  Its value, compared with a \\\chi^2_1\\.

- df:

  One, always.

- p.value:

  The upper tail of that \\\chi^2_1\\.

- estimate:

  The unrestricted estimate of the coefficient.

- null_value:

  The value it was tested against.

- parameter, coefficient:

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
refits nothing. They answer different questions and only the first is
printed.

`converged` is the inner optimizer's flag: whether a stopping rule
fired, a boolean about a threshold on a score whose size depends on the
model. `mode_error` is
[`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md),
how much of the penalized objective is still on the table at the point
the refit returned, in log-likelihood units against
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
– the rule statmodels7 0.81.0 established for exactly this question.

The flag is not merely a false negative: it is anti-correlated with the
quality of the point, because a stopping rule reads a change and a
distance reads a point. Measured on `y ~ x` with a Poisson response at
\\n = 200\\ over nineteen held values, all nineteen are at their mode
and the four the flag rejects are the four best-located of them. See
[`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md).

A statistic read off a refit that stopped short is the value at wherever
it stopped, and [`print()`](https://rdrr.io/r/base/print.html) says so –
reading `mode_error`, so that it says it where it is true.

The class exists rather than a bare list because everything a caller of
this toolkit receives is an object with declared properties. R's own
tests return an `htest`, which is a list with a class attribute and no
contract; nothing here can validate one, and a reader cannot ask it what
it carries.

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
S7::S7_inherits(statmod_test(fit, "mu", "x", 1), StatmodTest)
#> [1] TRUE
```
