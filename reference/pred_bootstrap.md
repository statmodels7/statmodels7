# Parametric Bootstrap Predictive Distribution

The parametric bootstrap predictive distribution of Harris (1989): `n`
responses are simulated from the fit with the random effects drawn
afresh, each response is refitted, and the family at each refit is
averaged over the effects. It does not rely on the normal approximation
of
[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
and it costs `n` refits.

## Usage

``` r
pred_bootstrap(n = 200L, refit = c("coefficients", "full"))
```

## Arguments

- n:

  The number of replicas, a single positive whole number. The default is
  `200`.

- refit:

  `"coefficients"` (the default), refitting each replica at the fit's
  hyperparameters, or `"full"`, choosing them again.

## Value

An S7 object of class `PredBootstrap`, inheriting from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md),
with the properties `label`, `n` (an integer) and `refit`.

## Details

With `refit = "coefficients"` each replica is refitted at the fit's
hyperparameters. With `refit = "full"` the hyperparameters are chosen
again on every replica, so their uncertainty enters the interval as
well, at the cost of a whole fit per replica. The bootstrap is rejected
for a model with a block that moves with its coefficients, for a
response that is not one number per observation, and, with
`refit = "full"`, for a cross-validated criterion. The replicas use the
random number generator of the session, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
reproduces them. See
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md)
and
[`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md).

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

## See also

[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## Examples

``` r
set.seed(2)
dd <- data.frame(x = runif(30))
dd$y <- 1 + dd$x + rnorm(30, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
predict(fit, what = "response", newdata = data.frame(x = 0.5),
        interval = prediction_interval(
          predictive = pred_bootstrap(n = 20, refit = "coefficients")))
#>        fit        se     lower    upper
#> 1 1.547902 0.5891418 0.3767825 2.695134
```
