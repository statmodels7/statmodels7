# Plug-in Predictive Distribution

The family at the estimates, averaged only over the effects that
`random` sets aside. The interval is the estimative one, which ignores
the error in the estimates and so covers less than the level in a small
sample. It costs one evaluation of the predictors and is the cheapest of
the three predictive objects.

## Usage

``` r
pred_plugin()
```

## Value

An S7 object of class `PredPlugin`, inheriting from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md),
with the one property `label`.

## Details

See
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## See also

[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
predict(fit, what = "response", newdata = data.frame(x = 0.5),
        interval = prediction_interval(predictive = pred_plugin()))
#>        fit        se     lower    upper
#> 1 1.529108 0.3337831 0.8749056 2.183311
```
