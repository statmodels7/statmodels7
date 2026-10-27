# Predict From a Fitted Model

Predicts from a fitted model: any one of the distribution's parameters,
any of its moments, every parameter or linear predictor at once, or a
new observation of the response, at the fitting data or at new data,
with standard errors and intervals on request.

## Usage

``` r
# S3 method for class 'StatmodFit'
predict(
  object,
  what = "parameter",
  newdata = NULL,
  se = FALSE,
  level = 0.95,
  random = "conditional",
  interval = confidence_interval(),
  vcov_type = bayesian_vcov(),
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what:

  What to predict: a parameter's name, optionally prefixed `"link:"`; a
  moment's name; `"parameter"` (the default) or `"link"`; or
  `"response"`, which requires `interval = prediction_interval()` and is
  the only value that this interval accepts. An unrecognized name
  signals an error listing what is available.

- newdata:

  A data frame, or `NULL` for the fitting data. It needs the variables
  of the equations that are read: where `what` names one parameter and
  `interval` is
  [`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
  the variables that only the other equations use may be absent. It does
  not need the response, except that a moment or a prediction interval
  of a binomial fitted with
  [`cbind()`](https://rdrr.io/r/base/cbind.html) reads the number of
  trials from it. For a model with a structural term, rows that all
  carry the response are read as a series of their own, and rows without
  it continue the fitted series. A data frame with no rows signals an
  error.

- se:

  `TRUE` to report the standard error and an interval as well. The
  default is `FALSE`. A moment has no standard error, and `se = TRUE`
  with a moment signals an error.
  [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
  sets it to `TRUE`, and
  [`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
  does not read it, because its result carries its own `se` column.

- level:

  The interval's level, `0.95` by default, a single number strictly
  between 0 and 1. It is read wherever an interval is returned: with
  `se = TRUE`, with
  [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
  and with
  [`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
  which gives its interval without `se`.

- random:

  How a
  [`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
  term is read: `"conditional"` (the default) with each group's own
  estimated effect, so a level that is absent from the fitting data
  signals an error; `"zero"` with every effect at zero, the typical
  group; `"marginal"` averaged over the prior that the fit estimated,
  the population average \\E_b\[h^{-1}(\eta + z^\top b)\]\\. A single
  string applies to every such term; a character vector named by the
  terms' keys sets the mode term by term, and the terms that it does not
  name stay conditional. See the section on new groups.

- interval:

  An interval object that selects the uncertainty that the interval
  describes.
  [`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md)
  (the default) gives the interval of the parameter, the predictor or
  the moment from the uncertainty of the estimates;
  [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
  gives the interval of the parameter of a new group;
  [`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
  gives the interval of a new observation of the response, and takes the
  predictive distribution as its argument
  (`predictive = pred_averaged()`,
  [`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
  or
  [`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)).
  A character string is rejected with an error. See
  [predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md)
  and
  [predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object that selects the variance convention of
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  (the default),
  [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  or
  [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md),
  which also accounts for the estimated smoothing parameters. It is read
  where `se` is `TRUE`, with
  [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
  and with
  [`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
  holding
  [`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
  and a value other than the default signals an error elsewhere. The
  unconditional variance differs most from the others in a band around a
  penalized term, since a smoothing parameter moves the fitted values.
  The conventions are compared at
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

- ...:

  Not used. Any argument signals an error that names it.

## Value

With `se = FALSE`, a numeric vector with one value per row of `newdata`
(of the fitting data where `newdata` is `NULL`) when `what` names one
quantity, and a named list of such vectors for `"parameter"` and
`"link"`.

With `se = TRUE`, a data frame with columns `fit`, `se`, `lower` and
`upper` in place of each vector, or, for `"parameter"` and `"link"`, a
named list of such data frames, one per parameter. `se` is `NA` for an
observation whose predictor reads a coefficient that has no variance.

With
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
a data frame with columns `fit` (the predictive median), `se` (the
predictive standard deviation, `NA` where it does not exist), `lower`
and `upper`, all `NA` for a row whose predictor reads a coefficient that
has no variance.

## Quantities that can be predicted

`what` takes

- a parameter's name:

  `"mu"`, `"sigma"` or `"alpha"`, whatever the family calls them. A
  parameter is available for every family, because the model fits it,
  and it exists even where a moment does not.

- a moment's name:

  `"mean"`, `"variance"`, `"std_dev"`, `"skewness"`, `"kurtosis"`. A
  moment is available where the family has one, and the result is `NaN`
  or `NA` where the moment does not exist. A Cauchy's mean is `NaN`,
  which is the correct answer for it. A moment has no standard error.
  For a binomial fitted with
  [`cbind()`](https://rdrr.io/r/base/cbind.html), a moment needs the
  number of trials of each row, so `newdata` has to carry the variables
  with which the response is written.

- `"parameter"`:

  every parameter at once, as a named list. This is the default.

- `"link"`:

  every linear predictor at once, before the inverse link.

- `"response"`:

  the response itself, with `interval = prediction_interval()`, which it
  requires: its median and the ends of an interval in which a new
  observation falls. See the section on intervals.

A parameter's name may be prefixed by `"link:"` to return its predictor
instead of its value, as `"link:sigma"`.

## The argument order departs from [`stats::predict()`](https://rdrr.io/r/stats/predict.html)

There the second argument is `newdata`; here it is `what`. A statmod fit
has several parameters and several moments, so choosing among them is
the ordinary variation and predicting on new data is the occasional one.

A data frame passed second is detected, and the error names the argument
order.

## New data

New data is passed through each term's blueprint, so a factor keeps the
levels and contrasts with which it was fitted, a spline its knots, and a
basis its reparametrization. The new frame supplies the values at which
the terms are evaluated and does not change how they are built. Offsets
are evaluated again at the new rows.

The new frame needs the variables of the equations that are read. Where
`what` names one parameter and `interval` is
[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
the other equations are not read, and the variables that only they use
may be absent from `newdata`. The response is not needed, except that a
binomial fitted with [`cbind()`](https://rdrr.io/r/base/cbind.html)
reads the number of trials from it for a moment or a prediction
interval.

## A score-driven term is predicted past the series

The contribution of such a term at one row is the state that a recursion
has reached, so new rows without the response continue the series
instead of being read on their own. Each row is placed by its own time
within its own group, and must come after every observed time of that
group. A row that falls inside the observed series is rejected, since
there the response is known and the filter has to be run, not continued.
A group that is absent from the fitting data is rejected as well, and a
term built without `time` cannot be continued, because the place of a
row in the series is then its position in the data frame.

New rows that carry the response are a series of their own: the
structural term is rebuilt on them with the fitted parameters, and they
are read as the fitting rows are, including the standard errors. If the
fitting data are passed as `newdata`, the result is that of
`predict(fit)`; this holds for a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term as well, whose prediction is the posterior-weighted predictor. A
frame that carries the response on some rows only is rejected. A
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term cannot be continued past the series, since what it reports at an
observed row is a posterior over states.

Beyond the data the score sits at its conditional mean of zero, which
the model's own definition guarantees, so the continuation is the
deterministic recursion and involves no simulation.

## A new group, and the typical one

A random-effect term predicts a group that is present in the fitting
data with that group's estimated effect. A group that is absent from the
fitting data has none, and what is put in its place depends on `random`.
`"zero"` sets every effect of the term to zero, at every row, which is
the prediction for the typical group (lme4's and glmmTMB's
`re.form = NA`). `"marginal"` averages over the prior instead,
\\E_b\[h^{-1}(\eta_0 + z^\top b)\]\\, which is the population average.
With an identity link the two agree for the mean; with a log link the
marginal mean is \\\exp(\eta_0 + \sigma_b^2/2)\\, and with a logit it is
pulled towards one half. A moment is averaged by the law of total
expectation and the variance adds the variance of the conditional mean,
so a random effect on the scale enters the marginal variance of the
response.

The average is a Gauss-Hermite product grid where every prior is
Gaussian and there are at most three coordinates in all; adaptive
quadrature against the prior's density
([`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html))
where the only effect set aside is one coordinate with a prior that is
not Gaussian; and 10000 draws from the priors otherwise, on a seed of
their own. A prior that is not Gaussian is integrated only where the
inverse link is bounded: a Student t has no moment generating function,
so under a log link the average does not exist, and such a request
signals an error. Under a link whose predictor has a restricted domain
(the square root and inverse links) the average does not exist either,
since a Gaussian effect leaves the domain with positive probability; it
is reported where that probability is below 1e-8 and is `NA`, with a
warning, elsewhere. A term that shares a covariance with others through
a label has no prior of its own, so `"marginal"` cannot average over it;
it can be read `"conditional"` or set aside at `"zero"`. The standard
error of a marginal parameter is the delta method on the fixed part,
conditional on the prior's scale.

A random effect written inside a term's subformula, such as
`nl(~ Asym / (1 + exp((xmid - age) / scal)), Asym ~ 1 + random(~ 1 | g))`
or `seg(t, psi ~ random(~ 1 | id))`, is set aside by the same modes. Its
key is the outer term's key, the parameter and the sub-term joined by
`"::"`, and an unnamed mode applies to it as well. The effect enters the
term's parameter, and the term is not linear in it, so the term is
evaluated at the effect's value rather than having its columns removed:
`"zero"` evaluates it with the effect at zero at every row, and the
delta method reads its Jacobian there. `"marginal"` averages the term's
contribution, carried through the inverse link, over the nodes that the
rule above builds for the priors set aside; a prior of another family
signals an error. A predictor (`"link"`) is the average of the
predictors at the nodes. The standard error of a marginal parameter is
the delta method of that average, whose gradient sums the design at each
node, and its interval is the delta method's on the link scale carried
through the link.
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
are not available for such an effect, and a standard error is not
available where an average is taken over effects of both kinds.

## Intervals and predictive distributions

The object given to `interval` selects the uncertainty that the interval
describes:
[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md)
the estimates alone,
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
the estimates and the effects of a new group, and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
a new observation of the response. The topic
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md)
compares the three and gives the formulas. A prediction interval holds a
predictive object that sets how the estimation uncertainty enters:
[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md)
(the default),
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
or
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md).
The topic
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md)
compares them and states what each costs.

[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
are rejected for a multivariate family and for a model with a structural
term. The structural term also rejects `random = "zero"` and
`random = "marginal"`: its level is a recursion read at the predictor,
so setting a random effect aside would move the recursion as well.

A forecast's standard error, with `se = TRUE`, is only the uncertainty
of the estimated parameters, carried by the delta method through the
continued recursion, whose derivative is exact. A forecast also carries
the uncertainty of the future scores, which the delta method does not
reach, so the standard error and the interval leave it out, and a
warning is issued for the equations that the continued filter reaches.

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

## See also

[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
for one parameter's fitted values,
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
for the matched diagnostic,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the variance from which the standard errors are computed.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)

# One parameter and one of the family's moments are predicted.
head(predict(fit, "mu"))
#> [1] 1.598797 1.804267 2.191116 2.837415 1.475789 2.818493
head(predict(fit, "variance"))
#> [1] 0.16546123 0.14306587 0.10879810 0.06885779 0.18051281 0.06978618

# The predictor of a parameter is returned instead of its value.
head(predict(fit, "link:sigma"))
#> [1] -0.8995092 -0.9722251 -1.1091307 -1.3378560 -0.8559768 -1.3311596

# For a Gaussian the mean is mu and the variance is sigma squared, and the
# moments of the fit return these values.
all.equal(predict(fit, "mean"), predict(fit, "mu"))
#> [1] TRUE
all.equal(predict(fit, "variance"), predict(fit, "sigma")^2)
#> [1] TRUE

# The result carries an interval as well.
head(predict(fit, "mu", se = TRUE))
#>        fit         se    lower    upper
#> 1 1.598797 0.07213967 1.457406 1.740188
#> 2 1.804267 0.05854267 1.689526 1.919009
#> 3 2.191116 0.04332413 2.106202 2.276029
#> 4 2.837415 0.06721771 2.705670 2.969159
#> 5 1.475789 0.08114458 1.316749 1.634830
#> 6 2.818493 0.06593854 2.689256 2.947730

# Every parameter is returned at once, on either scale.
str(predict(fit, "parameter"))
#> List of 2
#>  $ mu   : num [1:60] 1.6 1.8 2.19 2.84 1.48 ...
#>  $ sigma: num [1:60] 0.407 0.378 0.33 0.262 0.425 ...

# For a new group, the prediction is the typical one or the population
# average.
gg <- data.frame(g = factor(rep(letters[1:12], each = 10)),
                 x = rnorm(120))
gg$y <- rpois(120, exp(0.5 + 0.3 * gg$x + rnorm(12, sd = 0.6)[gg$g]))
fg <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(),
              gg)
new <- data.frame(g = "new", x = 0)
predict(fg, "mu", new, random = "zero")
#> [1] 1.593882
predict(fg, "mu", new, random = "marginal")
#> [1] 1.786905

# The calls return the interval of a new group's mean and of a new count
# from it.
predict(fg, what = "mu", newdata = new, random = "zero",
        interval = group_interval())
#>        fit        se     lower    upper
#> 1 1.593882 0.9711082 0.5944571 4.273577
predict(fg, what = "response", newdata = new, random = "zero",
        interval = prediction_interval(predictive = pred_plugin()))
#>   fit       se lower upper
#> 1   1 1.614653     0     6
predict(fg, what = "response", newdata = new, random = "zero",
        interval = prediction_interval())
#>   fit       se lower upper
#> 1   1 1.658936     0     6

# A name that the family does not have is rejected, and the message lists
# what is available.
try(predict(fit, "median"))
#> Error : 'median' is neither a parameter of this distribution nor a
#>   moment. The parameters are: mu, sigma.
#>   The moments are: mean, variance, std_dev, skewness, kurtosis.
#>   "parameter" and "link" give all of them at once.
```
