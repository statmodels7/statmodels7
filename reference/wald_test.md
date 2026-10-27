# Wald's Test of a Coefficient

Selects Wald's statistic for the hypothesis that one coefficient equals
a value \\b\\: \$\$W = \frac{(\hat\beta_j -
b)^2}{\widehat{\mathrm{Var}}(\hat\beta_j)},\$\$ compared with a
\\\chi^2_1\\ distribution. It reads the unrestricted fit and the
variance matrix of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and refits nothing, which is why it is the default statistic of every
function that takes a `test`.

## Usage

``` r
wald_test()
```

## Value

A `WaldTest` object, inheriting from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
with the properties `name` (`"wald"`) and `label` (`"Wald"`).

## Details

Choose it when the cost of a refit per coefficient matters, or when the
interval symmetric about the estimate is what is wanted. The likelihood
ratio, the score and the gradient statistics follow the shape of the
likelihood, and the Wald statistic does not; see
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison. The variance on which it depends is selected by the
argument `vcov_type` of the consuming function (see
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)).

The constructor has no arguments.

## See also

[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
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
wald_test()
#> <test> Wald
statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
             test = wald_test())
#> 
#>  Wald test on a statmod fit
#> 
#> coefficient:  mu:x
#> X-squared = 11.67,  df = 1,  p-value 0.0006347
#> alternative hypothesis: true x is not equal to 1.5
#> estimate: 2.008
confint(fit, test = wald_test())
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.043291 0.08770104  0.8714004
#> mu:x                     mu linpar           x  2.007856 0.14865515  1.7164973
#> sigma:(Intercept)     sigma linpar (Intercept) -1.032801 0.08006408 -1.1897235
#>                        upper
#> mu:(Intercept)     1.2151822
#> mu:x               2.2992148
#> sigma:(Intercept) -0.8758781
```
