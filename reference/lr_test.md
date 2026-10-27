# Likelihood-Ratio Test of a Coefficient

Selects the likelihood-ratio statistic for the hypothesis that one
coefficient equals a value \\b\\: \$\$LR = 2\\\ell_p(\hat\beta) -
\ell_p(\tilde\beta)\\,\$\$ where \\\ell_p\\ is the penalized
log-likelihood, \\\hat\beta\\ the estimates of the fit and
\\\tilde\beta\\ the estimates in which the coefficient is held at \\b\\.
It is compared with a \\\chi^2_1\\ distribution.

## Usage

``` r
lr_test()
```

## Value

An `LrTest` object, inheriting from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
with the properties `name` (`"lr"`) and `label` (`"Likelihood-ratio"`).

## Details

The statistic is a difference of the objective at two points, so a
reparametrization of the coefficient does not change it, and the
interval that inverts it follows the shape of the likelihood. It
requires one restricted refit per value \\b\\ (a few dozen per
interval), at the hyperparameters of the fit; where the outer criterion
estimated some of the coefficients, the unrestricted side is the joint
mode, obtained by one refit of the full model from the fitted
coefficients. See
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison with the other three.

The constructor has no arguments.

## See also

[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
and
[`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)
for the siblings,
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison of the four,
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
for the functions that take it.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
lr_test()
#> <test> Likelihood-ratio
statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
             test = lr_test())
#> 
#>  Likelihood-ratio test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 11.16,  df = 1,  p-value 0.0008379
#> alternative hypothesis: true x is not equal to 1.5
#> estimate: 2.008
confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE)
#>      parameter   term coefficient estimate        se    lower    upper
#> mu:x        mu linpar           x 2.007856 0.1486552 1.716674 2.299038
```
