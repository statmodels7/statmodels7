# Intervals of a Prediction

The `interval` argument of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
takes an interval object that selects the uncertainty that the interval
describes. There are three, and they differ in what the interval is
about.

- [`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md):

  the default, with `se = TRUE`: the interval of a parameter, a
  predictor or a moment, which carries only the uncertainty of the
  estimates. It costs one evaluation of the delta method.

- [`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md):

  the interval of the parameter of a new group, with `random = "zero"`
  or `"marginal"`. It adds the variation of the effects of a new group
  to the uncertainty of the estimates. It costs the delta method plus a
  quantile of the distribution of the predictor, by closed form under
  Gaussian priors and by quadrature or Monte Carlo otherwise.

- [`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md):

  the interval of a new observation of the response, with
  `what = "response"`. It adds the variation of the family itself. The
  way that the estimation uncertainty enters is set by its argument
  `predictive`, described in
  [predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).
  The cost is that of the predictive distribution chosen.

## Value

Not applicable: this topic describes the constructors. Each returns an
S7 object that inherits from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md).

## Confidence interval

With the default, an effect set aside by `random` gives the interval of
the typical group's value, and the value of a new group falls outside
that interval far more often than the level implies. The interval is
built on the scale of the equation and mapped back through the link.

## Group interval

A new group's parameter is \\\theta^\* = h^{-1}(\eta^\*)\\, with
\\\eta^\*\\ the typical group's predictor plus its estimation error plus
\\z^\top b\\, the effects of a new group drawn from the prior, in the
parameter's own equation, so a random effect on `sigma` widens the
interval for `sigma`. Under Gaussian priors the variance of \\\eta^\*\\
is \\se^2 + z^\top \Sigma_b z\\, the standard error that glmmTMB reports
at a new level on the link scale. The ends of the interval are quantiles
of \\\eta^\*\\ carried through the inverse link, which is monotone, so
they are the quantiles of \\\theta^\*\\ under any link; under a prior
that is not Gaussian they are the quantiles of the mixture of
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md).
Where the interval of \\\eta^\*\\ leaves the domain of the link (a
square root link and a negative predictor), the ends are `NA` at that
row, with a warning. The standard error is the standard deviation of
\\\theta^\*\\, not the delta method, and is `NA` where it is not known
to exist. The fit is the typical group's value under `random = "zero"`
and the population average under `"marginal"`. The interval implies
`se = TRUE`.

## Prediction interval

The family is averaged over the predictors of every parameter, which are
taken to be jointly Gaussian (correlated across equations where a label
ties them). The ends of the interval are quantiles of that average, and
the fit is its median, which every family has, so the family does not
need a location parameter. `se` is the standard deviation of the
average, and is `NA` where the family's variance does not exist. For a
discrete family the ends are values of the support, so the coverage is
at least the level rather than equal to it. Effects read `"conditional"`
enter with their own estimates, so this is also the interval of a new
observation of a group that is present in the fitting data.

## Restrictions

[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
are rejected for a multivariate family and for a model with a structural
term. The structural term also rejects `random = "zero"` and
`random = "marginal"`: its level is a recursion read at the predictor,
so setting a random effect aside would move the recursion as well.

Under a prior that is not Gaussian, both
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
read the quantiles of the prior: a Student t prior as a scale mixture of
Gaussians; a single univariate prior of another family, where it is the
only prior that is not Gaussian, on quantile nodes, which do not depend
on the seed; and any other combination of priors by Monte Carlo, which
does (see
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)).
A term that shares a label with an effect read with its own estimate, or
with one written inside a subformula, is rejected.

The level of every interval is the argument `level` of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## See also

[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
nd <- data.frame(x = c(0.2, 0.8))
predict(fit, what = "mu", newdata = nd, se = TRUE,
        interval = confidence_interval())
#>        fit         se    lower    upper
#> 1 1.283706 0.07987999 1.127144 1.440268
#> 2 1.774511 0.07581109 1.625924 1.923098
predict(fit, what = "response", newdata = nd,
        interval = prediction_interval())
#>        fit        se     lower    upper
#> 1 1.283706 0.3475094 0.5993255 1.968086
#> 2 1.774511 0.3465967 1.0918912 2.457131
```
