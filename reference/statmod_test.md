# Test One Coefficient Against One Value

Tests that one coefficient of a fit equals one value, by Wald's
statistic, the likelihood ratio, Rao's score or Terrell's gradient.

## Usage

``` r
statmod_test(
  fit,
  param,
  coefname,
  value = 0,
  test = wald_test(),
  vcov_type = bayesian_vcov()
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, coefname:

  Which coefficient, as in
  [`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md):
  the distribution parameter and the coefficient's name in its equation.

- value:

  The value under the null, a single number. The default is `0`, which
  is the value that a summary reports.

- test:

  A
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
  object:
  [`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md)
  (the default),
  [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
  [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
  or
  [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md).
  The first reads the unrestricted fit only; the other three cost one
  restricted refit. The four are compared at
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md).
  A character string is rejected with an error.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance that the Wald statistic reads, passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  (the default),
  [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  or
  [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md).
  The conventions are compared at
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

## Value

A
[`StatmodTest()`](https://statmodels7.github.io/statmodels7/reference/StatmodTest-class.md).

## Details

`summary(fit, test =)` reports one of the four statistics for the rows
of the coefficient tables that can be held (not a row under a kinked
penalty, not an aliased one, and none for a model with a structural
term), always against zero. This function tests one coefficient against
any value, for questions that a summary cannot answer, such as whether a
slope is one, whether an elasticity is unity, or whether a coefficient
matches a value fixed outside the data.

The statistics, the quantities that they read and their cost are
described at
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).
The last three are built on the restricted fit of
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
and inverting one gives the interval of
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

Reading one row costs one restricted refit, whereas `summary(test =)`
costs one per row. The function also reads the mode error of the refit,
which costs one Hessian of the restricted fit; the Rao score statistic
reads that matrix anyway, so that it adds no cost there. The summary and
the interval do not read the mode error.

## See also

[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the four statistics,
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the variance conventions,
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for the rows that can be held, against zero,
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
for the interval that a test inverts.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200))
dd$y <- rpois(200, exp(0.3 + dd$x))
fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
## the summary tests whether the slope is zero; this tests whether it is one
statmod_test(fit, param = "mu", coefname = "x", value = 1, test = lr_test())
#> 
#>  Likelihood-ratio test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 1.27,  df = 1,  p-value 0.2597
#> alternative hypothesis: true x is not equal to 1
#> estimate: 1.056
statmod_test(fit, param = "mu", coefname = "x", value = 0, test = lr_test())
#> 
#>  Likelihood-ratio test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 471.8,  df = 1,  p-value < 1e-16
#> alternative hypothesis: true x is not equal to 0
#> estimate: 1.056
```
