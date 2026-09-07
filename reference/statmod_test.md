# Test One Coefficient Against One Value

Tests that one coefficient of a fit equals one value, by Wald's
statistic, the likelihood ratio, Rao's score or Terrell's gradient.

## Usage

``` r
statmod_test(fit, param, coefname, value = 0, test = "wald", type = "bayesian")
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

  The value under the null, a single number. `0` by default, which is
  what a summary reports.

- test:

  One of `"wald"`, `"lr"`, `"score"`, `"gradient"`.

- type:

  Which variance the Wald statistic reads, passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A
[`StatmodTest()`](https://statmodels7.github.io/statmodels7/reference/StatmodTest-class.md).

## Details

`summary(fit, test =)` reports the same four statistics for every row of
the coefficient tables, and always against zero. This tests ONE
coefficient against ANY value, which is the question a summary cannot be
asked: whether a slope is one, whether an elasticity is unity, whether a
coefficient matches a value fixed outside the data.

What the statistics are, what they read and what they cost is at
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md);
the restricted fit the last three are built on is
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
and the interval obtained by inverting one is
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

Reading one row costs one restricted refit, where `summary(test =)`
costs one per row: measured on a Poisson fit of twenty-five covariates,
one statistic is 0.010 s against the summary's 0.220 s over twenty-six
rows. It also reads the refit's mode error, which is one Hessian – 63
per cent of the refit on that fit and 2.3 per cent on a penalized
smooth, and none at all for Rao's score, which reads the same matrix
anyway. The summary and the interval do not pay it.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for every row against zero,
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
for the interval a test inverts.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200))
dd$y <- rpois(200, exp(0.3 + dd$x))
fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
## the summary asks whether the slope is zero; this asks whether it is one
statmod_test(fit, "mu", "x", 1, "lr")
#> 
#>  Likelihood-ratio test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 1.27,  df = 1,  p-value 0.2597
#> alternative hypothesis: true x is not equal to 1
#> estimate: 1.056
statmod_test(fit, "mu", "x", 0, "lr")
#> 
#>  Likelihood-ratio test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 471.8,  df = 1,  p-value < 1e-16
#> alternative hypothesis: true x is not equal to 0
#> estimate: 1.056
```
