# The Confidence Interval of a Prediction

Selects the interval of a parameter, a predictor or a moment that
carries only the uncertainty of the estimates. It is the default
`interval` of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
and is reported with `se = TRUE`. Choose it for the typical group's
value; for the value of a new group use
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
and for a new response use
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md).

## Usage

``` r
confidence_interval()
```

## Value

An S7 object of class `ConfidenceInterval`, inheriting from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md),
with the one property `label`.

## Details

The standard error is the delta method on the scale of the equation, and
the interval is carried through the link. For an effect set aside by
`random` the interval is that of the typical group's value. See
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md).

## See also

[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
predict(fit, what = "mu", newdata = data.frame(x = 0.5), se = TRUE,
        interval = confidence_interval())
#>        fit         se    lower    upper
#> 1 1.529108 0.05284833 1.425528 1.632689
```
