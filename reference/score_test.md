# Rao's Score Test of a Coefficient

Selects Rao's score statistic for the hypothesis that one coefficient
equals a value \\b\\: \$\$S =
U_j(\tilde\beta)^2\\\[K(\tilde\beta)^{-1}\]\_{jj},\$\$ where
\\\tilde\beta\\ are the estimates in which the coefficient is held at
\\b\\, \\U_j\\ is the penalized score of the coefficient and \\K\\ the
penalized information. It is compared with a \\\chi^2_1\\ distribution.

## Usage

``` r
score_test()
```

## Value

A `ScoreTest` object, inheriting from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
with the properties `name` (`"score"`) and `label` (`"Rao score"`).

## Details

The statistic is evaluated at the restricted estimates, so the
unrestricted fit is not read. It requires one restricted refit per value
\\b\\, at the hyperparameters of the fit, and the inverse of the
penalized information at the restricted point, so each evaluation also
costs one Hessian. The general form \\U'K^{-1}U\\ reduces to the
expression above because every other component of the score vanishes at
the restricted maximum. See
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison with the other three.

The constructor has no arguments.

## See also

[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md)
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
score_test()
#> <test> Rao score
statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
             test = score_test())
#> 
#>  Rao score test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 10.41,  df = 1,  p-value 0.001252
#> alternative hypothesis: true x is not equal to 1.5
#> estimate: 2.008
confint(fit, parm = "mu:x", test = score_test(), readable = FALSE)
#>      parameter   term coefficient estimate        se    lower    upper
#> mu:x        mu linpar           x 2.007856 0.1486552 1.712996 2.302716
```
