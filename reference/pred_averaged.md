# Predictive Distribution Averaged Over the Estimates

The default `predictive` of
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md):
the family averaged over the predictors at the fit's estimates, with the
covariance of the estimates plus the covariance that the set-aside
effects add. It is the normal approximation to the parametric bootstrap
predictive distribution (Harris, 1989), and it covers a new observation
at close to the level in a small sample, where
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
covers less. It costs one evaluation of the covariance of the estimates,
so it is the choice when
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
is too expensive.

## Usage

``` r
pred_averaged()
```

## Value

An S7 object of class `PredAveraged`, inheriting from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md),
with the one property `label`.

## Details

It is the only predictive that reads the argument `type` of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).
See
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

## See also

[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md),
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
        interval = prediction_interval(predictive = pred_averaged()))
#>        fit        se     lower   upper
#> 1 1.529108 0.3423082 0.8547573 2.20346
```
