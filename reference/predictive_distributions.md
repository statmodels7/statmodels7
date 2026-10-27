# Predictive Distributions of a New Response

The argument `predictive` of
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
takes an object that selects how the estimation uncertainty enters the
distribution of a new response. There are three.

- [`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md):

  the default: the family averaged over the approximate sampling
  distribution of the estimates. It costs one evaluation of the
  covariance of the estimates.

- [`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md):

  the family at the estimates, ignoring their error. It costs one
  evaluation of the predictors.

- [`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md):

  the parametric bootstrap predictive distribution, by simulation and
  refit of `n` responses. It costs `n` refits, each as expensive as a
  fit with the hyperparameters held (`refit = "coefficients"`) or as a
  whole fit (`refit = "full"`).

## Value

Not applicable: this topic describes the constructors. Each returns an
S7 object that inherits from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md).

## Details

A new observation varies for three reasons: the family itself, the
effects of a new group where `random` sets them aside, and the error in
the estimates. The first two are random in the model and always enter.
The third enters according to the predictive object, for every parameter
alike, whether the family has one parameter (a Poisson) or several, and
whether or not a parameter has covariates.

## Averaged

The family is averaged over the predictors at the fit's estimates with
the covariance of the estimates from
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
plus the covariance that the effects add. This is the normal
approximation to the parametric bootstrap predictive distribution of
Harris (1989): the estimates are averaged over their approximate
sampling distribution, centered at the fit. For a Gaussian with a
constant \\\sigma = e^{\gamma_0}\\, \\\hat\gamma_0\\ of variance \\v\\,
the variance of a new observation is \\\mathrm{se}\_0^2 + z^\top\Sigma_b
z + \hat\sigma^2 e^{2v}\\, where the last term is the average of
\\\hat\sigma^{\*2}\\ over that distribution. The argument `type` of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the variance type of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
is read only with this predictive.

## Plug-in

The family at the estimates, averaged only over the effects: the
estimative interval, which ignores the error in the estimates and so
covers less than the level in a small sample.

## Bootstrap

The parametric bootstrap predictive distribution itself: `n` responses
are simulated from the fit with the random effects drawn afresh, each
response is refitted, and the family at each refit is averaged over the
effects. `refit = "coefficients"` refits at the fit's hyperparameters;
`"full"` chooses them again on every replica, so their uncertainty
enters too, at the cost of a whole fit each. The bootstrap is rejected
for a model with a block that moves with its coefficients
([`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)),
for a response that is not one number per observation, and, with
`refit = "full"`, for a cross-validated criterion. A replica whose refit
fails is left out, with a warning that reports how many were. The seed
is the caller's: [`set.seed()`](https://rdrr.io/r/base/Random.html)
before the call fixes the replicas. See
[`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md).

## Hyperparameters

The interval for a new group
([`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md))
is conditional on the covariance that the fit estimated, and so is the
interval for a new observation, so the uncertainty of the
hyperparameters is not propagated. The exception is
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
with `refit = "full"`, where the hyperparameters are chosen again on
every replica.

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

## See also

[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md),
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
nd <- data.frame(x = 0.5)
for (p in list(pred_averaged(), pred_plugin())) {
  print(predict(fit, what = "response", newdata = nd,
                interval = prediction_interval(predictive = p)))
}
#>        fit        se     lower   upper
#> 1 1.529108 0.3423082 0.8547573 2.20346
#>        fit        se     lower    upper
#> 1 1.529108 0.3337831 0.8749056 2.183311
```
