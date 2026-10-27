# The Interval of a New Observation

Selects the interval in which a new observation of the response falls.
It requires `what = "response"` in
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
and does not read `se`, because its result carries its own `se` column.
It adds the variation of the family to that of the estimates and of the
effects set aside; the argument `predictive` sets how the estimation
uncertainty enters.

## Usage

``` r
prediction_interval(predictive = pred_averaged())
```

## Arguments

- predictive:

  A
  [`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md)
  object:
  [`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md)
  (the default) averages the family over the approximate sampling
  distribution of the estimates,
  [`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
  uses the estimates alone, and
  [`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
  simulates and refits. See
  [predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## Value

An S7 object of class `PredictionInterval`, inheriting from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md),
with the properties `label` and `predictive`.

## Details

The family is averaged over the predictors of every parameter. See
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md)
for the interval and
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md)
for the three values of `predictive`.

## See also

[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md),
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
nd <- data.frame(x = 0.5)
predict(fit, what = "response", newdata = nd,
        interval = prediction_interval())
#>        fit        se     lower   upper
#> 1 1.529108 0.3423082 0.8547573 2.20346
predict(fit, what = "response", newdata = nd,
        interval = prediction_interval(predictive = pred_plugin()))
#>        fit        se     lower    upper
#> 1 1.529108 0.3337831 0.8749056 2.183311
```
