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
  random = "conditional",
  interval = c("confidence", "group", "prediction"),
  predictive = c("averaged", "plugin", "bootstrap"),
  n_boot = 200L,
  boot_refit = c("coefficients", "full"),
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

- random:

  How a
  [`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
  term is read: `"conditional"` (the default) with each group's own
  estimated effect, so a level the fit never saw signals an error;
  `"zero"` with every effect at zero, the typical group; `"marginal"`
  averaged over the prior the fit estimated, the population average
  \\E_b\[h^{-1}(\eta + z^\top b)\]\\. A single string applies to every
  such term; a character vector named by the terms' keys chooses term by
  term, the rest staying conditional. See the section on new groups.

- interval:

  `"confidence"` (the default), `"group"` or `"prediction"`. See the
  section on intervals.

- predictive:

  With `interval = "prediction"`, what the family is averaged over:
  `"averaged"` (the default), `"plugin"` or `"bootstrap"`. See the
  section on the predictive distribution.

- n_boot:

  The number of bootstrap replicas, `200` by default. Read only where
  `predictive = "bootstrap"`.

- boot_refit:

  `"coefficients"` (the default), refitting each replica at the fit's
  hyperparameters, or `"full"`, choosing them again. Read only where
  `predictive = "bootstrap"`.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  where `se` is `TRUE`. That is where `type` chooses between the
  Bayesian variance, the frequentist one and the unconditional one. A
  band around a penalized term is where the last of the three differs
  most from the others, since it is the fitted values that a smoothing
  parameter moves: measured on a univariate smooth at \\n = 200\\,
  `type = "unconditional"` widens this interval by 1.1 per cent on
  average and 7.6 per cent at its widest point.

## Value

With `se = FALSE`, a numeric vector of `nrow(newdata)` values when
`what` names one quantity, and a named list of such vectors for
`"parameter"` and `"link"`.

With `se = TRUE`, a data frame with columns `fit`, `se`, `lower` and
`upper` in place of each vector. `se` is `NA` for an observation whose
predictor reads a coefficient that has no variance, which is the truth
about such a fit, and no gap in the arithmetic.

With `interval = "prediction"`, a data frame with columns `fit` (the
predictive median), `se` (the predictive standard deviation), `lower`
and `upper`.

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

- `"response"`:

  the response itself, with `interval = "prediction"`: its median and
  the ends of an interval that a new observation falls in. See the
  section on intervals.

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
reached, so new rows without the response continue the series instead of
being read on their own. Each row is placed by its own time within its
own group, and must come after every observed time of that group; a row
falling inside the observed series is rejected, since there the response
is known and the filter must be run, never continued.

New rows that carry the response are a series of their own: the
structural term is rebuilt on them with the fitted parameters, and they
are read as the fitting rows are, standard errors included. With the
fitting data as `newdata` the result is `predict(fit)`'s; this holds for
a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term as well, whose prediction is the posterior-weighted predictor.

Beyond the data the score sits at its conditional mean of zero, which
the model's own definition guarantees, so the continuation is the
deterministic recursion and involves no simulation.

## A new group, and the typical one

A random-effect term predicts a group the fit saw with that group's
estimated effect. A group it never saw has none, and `random` says what
to put in its place. `"zero"` sets every effect of the term to zero, at
every row, which is the prediction for the typical group (lme4's and
glmmTMB's `re.form = NA`). `"marginal"` averages over the prior instead,
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
a label has no prior of its own and is read at `"zero"` only. The
standard error of a marginal parameter is the delta method on the fixed
part, conditional on the prior's scale.

A random effect written inside a term's subformula, such as
`nl(~ Asym / (1 + exp((xmid - age) / scal)), Asym ~ 1 + random(~ 1 | g))`
or `seg(t, psi ~ random(~ 1 | id))`, is set aside by the same modes. Its
key is the outer term's key, the parameter and the sub-term joined by
`"::"`, and an unnamed mode applies to it as well. The effect enters the
term's parameter, and the term is not linear in it, so the term is
evaluated at the effect's value rather than having its columns removed:
`"zero"` evaluates it with the effect at zero at every row, and the
delta method reads its Jacobian there. `"marginal"` averages the term's
contribution, carried through the inverse link, over the nodes of the
product grid of a Gaussian prior; a prior of another family signals an
error. A predictor (`"link"`) is the average of the predictors at the
nodes. The standard error of a marginal parameter is the delta method of
that average, whose gradient sums the design at each node, and its
interval is the delta method's on the link scale carried through the
link. `interval = "group"` and `"prediction"` are not available for such
an effect, and neither is a standard error where effects of both kinds
are averaged over.

## Three intervals

`interval` says what an interval is for.

- `"confidence"`:

  the default, with `se = TRUE`: the uncertainty of the estimates alone.
  For a term set aside by `random` it is the interval of the typical
  group's value, which a new group's value falls outside of far more
  often than the level says.

- `"group"`:

  the interval of a new group's parameter, with `random = "zero"` or
  `"marginal"`. A new group's parameter is \\\theta^\* =
  h^{-1}(\eta^\*)\\, with \\\eta^\*\\ the typical group's predictor plus
  its estimation error plus \\z^\top b\\, the effects of a new group
  drawn from the prior, in the parameter's own equation, so a random
  effect on `sigma` widens the interval for `sigma`. Under Gaussian
  priors the variance of \\\eta^\*\\ is \\se^2 + z^\top \Sigma_b z\\,
  the standard error glmmTMB reports at a new level on the link scale.
  The interval's ends are quantiles of \\\eta^\*\\ carried through the
  inverse link, which is monotone, so they are the quantiles of
  \\\theta^\*\\ under any link; under a prior that is not Gaussian they
  are the quantiles of the mixture of
  [`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md).
  The standard error is the standard deviation of \\\theta^\*\\, not the
  delta method, and is `NA` where it is not known to exist (see
  [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)).
  The fit is the typical group's value under `random = "zero"` and the
  population average under `"marginal"`. Implies `se = TRUE`.

- `"prediction"`:

  an interval for a new observation of the response, with
  `what = "response"`. The family is averaged over the predictors of
  every parameter, taken jointly Gaussian – correlated across equations
  where a label ties them – and the interval's ends are quantiles of
  that average and the fit is its median, which every family has, so no
  location parameter is needed; `se` is the standard deviation of the
  average, `NA` where the family's variance does not exist. For a
  discrete family the ends are values of the support, so the coverage is
  at least the level rather than equal to it. Effects read
  `"conditional"` enter with their own estimates, so this is also the
  interval of a new observation of a group the fit saw. What the average
  is over is `predictive`, below.

## The predictive distribution

A new observation varies for three reasons: the family itself, the
effects of a new group where `random` sets them aside, and the error in
the estimates. The first two are random in the model and always enter.
The third is what `predictive` decides, for every parameter alike,
whether the family has one parameter (a Poisson) or several, modelled or
not.

- `"averaged"`:

  the default. The family is averaged over the predictors at the fit's
  estimates with the covariance of the estimates from
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  plus the one the effects add. This is the normal approximation to the
  parametric bootstrap predictive distribution of Harris (1989): the
  estimates are averaged over their approximate sampling distribution,
  centred at the fit. For a Gaussian with a constant \\\sigma =
  e^{\gamma_0}\\, \\\hat\gamma_0\\ of variance \\v\\, the variance of a
  new observation is \\\mathrm{se}\_0^2 + z^\top\Sigma_b z +
  \hat\sigma^2 e^{2v}\\, the last term being the average of
  \\\hat\sigma^{\*2}\\ over that distribution.

- `"plugin"`:

  the family at the estimates, averaged over the effects alone: the
  estimative interval, which ignores the error in the estimates and so
  covers less than the level in a small sample.

- `"bootstrap"`:

  the parametric bootstrap predictive distribution itself: `n_boot`
  responses simulated from the fit, with the random effects drawn
  afresh, each refitted, and the family at each refit averaged over the
  effects. `boot_refit = "coefficients"` refits at the fit's
  hyperparameters; `"full"` chooses them again on every replica, so
  their uncertainty enters too, at the cost of a whole fit each. See
  [`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md).

Measured at 95 per cent over 400 samples, the coverage of a new
observation is 0.887 (plug-in) against 0.925 (averaged) for a Gaussian
regression at \\n = 12\\ (over a further 200 samples, 0.900, 0.950 and
0.945 with the bootstrap), and 0.935 against 0.948 for a Gamma with its
dispersion modelled at \\n = 25\\. Both intervals for a new group are
conditional on the covariance the fit estimated: its uncertainty is not
propagated. Measured at 95 per cent over 200 fits of a random intercept
at eight observations a group, a new group's predictor is covered 0.898
of the time over 10 groups and 0.943 over 40 by `"group"`, and a new
observation 0.939 and 0.951 by `"prediction"`, where the confidence
interval of the typical group covers 0.45 and 0.26; on a Poisson
response the same read 0.909 and 0.937 and, the support being discrete,
0.975 and 0.976. Under a prior that is not Gaussian both
`interval = "group"` and `interval = "prediction"` read its quantiles, a
Student t prior as a scale mixture of Gaussians and any other prior by
Monte Carlo (see
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)).
A term that shares a label with an effect read with its own estimate, or
with one written inside a subformula, is rejected.

A forecast's standard error, with `se = TRUE`, is the uncertainty of the
estimated parameters alone, carried by the delta method through the
continued recursion, whose derivative is exact. A forecast also carries
the uncertainty of the future scores, which the delta method does not
reach, so the standard error and the interval leave it out, and a
warning says so.

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

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
#> [1] 1.598797 1.804267 2.191116 2.837415 1.475789 2.818493
head(predict(fit, "variance"))
#> [1] 0.16546123 0.14306587 0.10879810 0.06885779 0.18051281 0.06978618

# A parameter's predictor instead of its value.
head(predict(fit, "link:sigma"))
#> [1] -0.8995092 -0.9722251 -1.1091307 -1.3378560 -0.8559768 -1.3311596

# For a Gaussian the mean is mu and the variance is sigma squared, which
# is what the moments come to.
all.equal(predict(fit, "mean"), predict(fit, "mu"))
#> [1] TRUE
all.equal(predict(fit, "variance"), predict(fit, "sigma")^2)
#> [1] TRUE

# With an interval.
head(predict(fit, "mu", se = TRUE))
#>        fit         se    lower    upper
#> 1 1.598797 0.07204781 1.457586 1.740008
#> 2 1.804267 0.05848939 1.689630 1.918904
#> 3 2.191116 0.04332275 2.106204 2.276027
#> 4 2.837415 0.06713927 2.705824 2.969005
#> 5 1.475789 0.08102934 1.316975 1.634604
#> 6 2.818493 0.06586368 2.689403 2.947584

# Every parameter at once, on either scale.
str(predict(fit, "parameter"))
#> List of 2
#>  $ mu   : num [1:60] 1.6 1.8 2.19 2.84 1.48 ...
#>  $ sigma: num [1:60] 0.407 0.378 0.33 0.262 0.425 ...

# A new group: the typical one, and the population average.
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

# The interval of a new group's mean, and of a new count from it.
predict(fg, "mu", new, random = "zero", interval = "group")
#>        fit        se     lower    upper
#> 1 1.593882 0.9711082 0.5944571 4.273577
predict(fg, "response", new, random = "zero", interval = "prediction",
        predictive = "plugin")
#>   fit       se lower upper
#> 1   1 1.614653     0     6
predict(fg, "response", new, random = "zero", interval = "prediction")
#>   fit       se lower upper
#> 1   1 1.658936     0     6

# A name the family does not have is refused, and the message says what
# is available.
try(predict(fit, "median"))
#> Error : 'median' is neither a parameter of this distribution nor a
#>   moment. The parameters are: mu, sigma.
#>   The moments are: mean, variance, std_dev, skewness, kurtosis.
#>   "parameter" and "link" give all of them at once.
```
