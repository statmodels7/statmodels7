# Predict From a Fitted Model

Any one of the distribution's parameters, any of its moments, or all of
the parameters or linear predictors at once, at the fitting data or at
new data.

## Usage

``` r
# S3 method for class 'StatmodFit'
predict(object, what = "parameter", newdata = NULL, ...)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what:

  What to predict: a parameter's name, a moment's name, `"parameter"` or
  `"link"`.

- newdata:

  A data frame, or `NULL` for the fitting data.

- ...:

  Unused.

## Value

A numeric vector when `what` names one quantity, and a named list of
vectors for `"parameter"` and `"link"`.

## Details

**What can be asked for.** `what` takes

- a parameter's name:

  `"mu"`, `"sigma"`, `"alpha"` – whatever the family calls them. Always
  available, whatever the family: a parameter is what the model fits,
  and it exists even where a moment does not.

- a moment's name:

  `"mean"`, `"variance"`, `"std_dev"`, `"skewness"`, `"kurtosis"`.
  Available where the family has one, and answering `NaN` or `NA` where
  it does not exist – a Cauchy's mean is `NaN`, which is the honest
  answer and not a failure.

- `"parameter"`:

  every parameter at once, as a named list. The default.

- `"link"`:

  every linear predictor at once, before the inverse link.

A parameter's name may be prefixed by `"link:"` to ask for its predictor
instead of its value, as `"link:sigma"`.

**The argument order departs from
[`predict`](https://rdrr.io/r/stats/predict.html)**, where the second
argument is `newdata`. Here it is `what`, because a statmod fit has
several parameters and several moments and choosing among them is the
ordinary variation, while predicting on new data is the occasional one.
Passing a data frame second is caught and named rather than failing
somewhere inside.

**New data** goes through the terms' blueprints, so a factor keeps the
levels and the contrasts it was fitted with rather than being rebuilt
from whatever the new frame happens to contain.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`fitted.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
head(predict(fit, "mu"))
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('statmodels7::StatmodFit', 'S7_object')"
head(predict(fit, "variance"))
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('statmodels7::StatmodFit', 'S7_object')"
head(predict(fit, "link:sigma"))
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('statmodels7::StatmodFit', 'S7_object')"
```
