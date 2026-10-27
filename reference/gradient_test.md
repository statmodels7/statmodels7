# Terrell's Gradient Test of a Coefficient

Selects Terrell's gradient statistic for the hypothesis that one
coefficient equals a value \\b\\: \$\$T =
U_j(\tilde\beta)\\(\hat\beta_j - b),\$\$ where \\\hat\beta\\ are the
estimates of the fit, \\\tilde\beta\\ the estimates in which the
coefficient is held at \\b\\ and \\U_j\\ the penalized score of the
coefficient. It is compared with a \\\chi^2_1\\ distribution.

## Usage

``` r
gradient_test()
```

## Value

A `GradientTest` object, inheriting from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
with the properties `name` (`"gradient"`) and `label`
(`"Terrell gradient"`).

## Details

The statistic is the product of the score and the distance between the
two fits. It requires one restricted refit per value \\b\\, at the
hyperparameters of the fit. Unlike the score statistic it reads no
Hessian and inverts no matrix, so a refit is its whole cost. See
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison with the other three.

The constructor has no arguments.

## See also

[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md)
and
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
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
gradient_test()
#> <test> Terrell gradient
statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
             test = gradient_test())
#> 
#>  Terrell gradient test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 10.41,  df = 1,  p-value 0.001252
#> alternative hypothesis: true x is not equal to 1.5
#> estimate: 2.008
confint(fit, parm = "mu:x", test = gradient_test(), readable = FALSE)
#>      parameter   term coefficient estimate        se    lower    upper
#> mu:x        mu linpar           x 2.007856 0.1486552 1.712996 2.302716
```
