# Predict From a Fitted Model

Any one of the distribution's parameters, any of its moments, or all of
the parameters or linear predictors at once, at the fitting data or at
new data.

## Usage

``` r
# S3 method for class 'StatmodFit'
predict(
  object,
  what = "parameter",
  newdata = NULL,
  se = FALSE,
  level = 0.95,
  ...
)
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

- se:

  Whether to report the prediction's uncertainty as well.

- level:

  The interval's level, where `se` is `TRUE`.

- ...:

  Passed to
  [`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  where `se` is `TRUE`, which is where `type` chooses between the
  bayesian variance and the frequentist one.

## Value

A numeric vector when `what` names one quantity, and a named list of
vectors for `"parameter"` and `"link"`. With `se = TRUE`, a data frame
of `fit`, `se`, `lower` and `upper` in place of each vector.

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

**A model carrying a score-driven term is predicted past the series.**
Such a term's contribution at one row is the state a recursion has
reached, so new rows continue the series rather than being read on their
own: each row is placed by its own time within its own group, and must
come after every observed time of that group. Beyond the data the score
sits at its conditional mean of zero, which the model's own definition
guarantees, so the continuation is the deterministic recursion and
involves no simulation. A forecast reports no standard error –
`se = TRUE` gives the uncertainty of the parameters, while a forecast
carries the uncertainty of the future scores as well, which is the
larger part and is no delta method.

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
#> [1] 1.598880 1.804324 2.191124 2.837341 1.475888 2.818422
head(predict(fit, "variance"))
#> [1] 0.16066196 0.13864759 0.10505462 0.06608504 0.17548015 0.06698799
head(predict(fit, "link:sigma"))
#> [1] -0.9142264 -0.9879100 -1.1266374 -1.3584065 -0.8701147 -1.3516210
head(predict(fit, "mu", se = TRUE))
#>        fit         se    lower    upper
#> 1 1.598880 0.07098561 1.459751 1.738009
#> 2 1.804324 0.05761990 1.691391 1.917257
#> 3 2.191124 0.04258405 2.107660 2.274587
#> 4 2.837341 0.06585185 2.708274 2.966408
#> 5 1.475888 0.07983104 1.319422 1.632354
#> 6 2.818422 0.06459871 2.691811 2.945033
```
