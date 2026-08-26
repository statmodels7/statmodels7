# Predict From a Fitted Model

Predicts from a fitted model: any one of the distribution's parameters,
any of its moments, or every parameter or linear predictor at once, at
the fitting data or at new data, with standard errors and intervals on
request.

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
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what:

  What to predict: a parameter's name, optionally prefixed `"link:"`; a
  moment's name; `"parameter"` (the default) or `"link"`. An
  unrecognized name signals an error listing what is available.

- newdata:

  A data frame, or `NULL` for the fitting data. Needs the covariates the
  model names but not the response.

- se:

  `TRUE` to report the standard error and an interval as well. `FALSE`
  by default.

- level:

  The interval's level, `0.95` by default. Read only where `se` is
  `TRUE`.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  where `se` is `TRUE`. That is where `type` chooses between the
  Bayesian variance and the frequentist one.

## Value

With `se = FALSE`, a numeric vector of `nrow(newdata)` values when
`what` names one quantity, and a named list of such vectors for
`"parameter"` and `"link"`.

With `se = TRUE`, a data frame with columns `fit`, `se`, `lower` and
`upper` in place of each vector. `se` is `NA` for an observation whose
predictor reads a coefficient that has no variance, which is the truth
about such a fit, and no gap in the arithmetic.

## What can be asked for

`what` takes

- a parameter's name:

  `"mu"`, `"sigma"`, `"alpha"` – whatever the family calls them. Always
  available, whatever the family: a parameter is what the model fits,
  and it exists even where a moment does not.

- a moment's name:

  `"mean"`, `"variance"`, `"std_dev"`, `"skewness"`, `"kurtosis"`.
  Available where the family has one, and answering `NaN` or `NA` where
  it does not exist. A Cauchy's mean is `NaN`, which is the correct
  answer for it.

- `"parameter"`:

  every parameter at once, as a named list. The default.

- `"link"`:

  every linear predictor at once, before the inverse link.

A parameter's name may be prefixed by `"link:"` to ask for its predictor
instead of its value, as `"link:sigma"`.

## The argument order departs from [`stats::predict()`](https://rdrr.io/r/stats/predict.html)

There the second argument is `newdata`; here it is `what`. A statmod fit
has several parameters and several moments, so choosing among them is
the ordinary variation and predicting on new data is the occasional one.

Passing a data frame second is caught and named, instead of failing
somewhere inside.

## New data

Goes through each term's blueprint, so a factor keeps the levels and
contrasts it was fitted with, a spline its knots, and a basis its
reparametrization. Nothing is rebuilt from whatever the new frame
happens to contain.

## A score-driven term is predicted past the series

Such a term's contribution at one row is the state a recursion has
reached, so new rows continue the series instead of being read on their
own. Each row is placed by its own time within its own group, and must
come after every observed time of that group; a row falling inside the
observed series is refused, since there the response is known and the
filter must be run, never continued.

Beyond the data the score sits at its conditional mean of zero, which
the model's own definition guarantees, so the continuation is the
deterministic recursion and involves no simulation.

A forecast reports **no standard error**. `se = TRUE` gives the
uncertainty of the parameters, while a forecast carries the uncertainty
of the future scores as well, which is the larger part and is no delta
method. Reporting the smaller half under the name of the whole would
mislead.

## See also

[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
for one parameter's fitted values,
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
for the matched diagnostic,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the variance the standard errors come from.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)

# One parameter, and one of the family's moments.
head(predict(fit, "mu"))
#> [1] 1.598880 1.804324 2.191124 2.837341 1.475888 2.818422
head(predict(fit, "variance"))
#> [1] 0.16066196 0.13864759 0.10505462 0.06608504 0.17548015 0.06698799

# A parameter's predictor instead of its value.
head(predict(fit, "link:sigma"))
#> [1] -0.9142264 -0.9879100 -1.1266374 -1.3584065 -0.8701147 -1.3516210

# For a Gaussian the mean is mu and the variance is sigma squared, which
# is what the moments come to.
all.equal(predict(fit, "mean"), predict(fit, "mu"))
#> [1] TRUE
all.equal(predict(fit, "variance"), predict(fit, "sigma")^2)
#> [1] TRUE

# With an interval.
head(predict(fit, "mu", se = TRUE))
#>        fit         se    lower    upper
#> 1 1.598880 0.07098561 1.459751 1.738009
#> 2 1.804324 0.05761990 1.691391 1.917257
#> 3 2.191124 0.04258405 2.107660 2.274587
#> 4 2.837341 0.06585185 2.708274 2.966408
#> 5 1.475888 0.07983104 1.319422 1.632354
#> 6 2.818422 0.06459871 2.691811 2.945033

# Every parameter at once, on either scale.
str(predict(fit, "parameter"))
#> List of 2
#>  $ mu   : num [1:60] 1.6 1.8 2.19 2.84 1.48 ...
#>  $ sigma: num [1:60] 0.401 0.372 0.324 0.257 0.419 ...

# A name the family does not have is refused, and the message says what
# is available.
try(predict(fit, "median"))
#> Error : 'median' is neither a parameter of this distribution nor a
#>   moment. The parameters are: mu, sigma.
#>   The moments are: mean, variance, std_dev, skewness, kurtosis.
#>   "parameter" and "link" give all of them at once.
```
