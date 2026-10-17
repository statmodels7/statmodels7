# The Parametric Bootstrap of a Prediction Interval

Returns the components of the predictive mixture that
[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md)
reads for `predictive = "bootstrap"`. There is one component per
replica, which is the fit refitted on a response simulated from the
fitted model and read at the rows predicted with no estimation
uncertainty of its own.

## Usage

``` r
predictive_bootstrap(object, spec, design, aside, n_boot, refit)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  The specification and design at the rows predicted, in which the
  columns of the effects set aside are zero.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that are not `"conditional"`.

- n_boot:

  The number of replicas.

- refit:

  `"coefficients"` or `"full"`.

## Value

A list of components, each a list with `eta`, `C` and `weight`. The list
has the attribute `infinite_variance`, which is `TRUE` where a prior has
no variance. The weights of the components add up to one.

## Details

A replica draws the response from the family at the fitted model. The
effects of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term that `random` sets aside are drawn afresh from the prior that the
fit estimated, so the groups of the data are new groups in every
replica, as they are in the sampling distribution of the
population-level estimates on which a new group's prediction rests. The
effects of a term read `"conditional"` keep their estimates: the
prediction is about those groups, and their effects are the quantity
predicted rather than a source of variation. Every other coefficient
keeps its estimate. The replica is then refitted in one of two ways:

- `"coefficients"`:

  every coefficient is re-estimated at the hyperparameters that the fit
  reached. A coefficient that the fit estimated on its criterion (a
  dispersion under
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md))
  is estimated there again while the hyperparameters are held, so the
  replica's estimator is the fit's.

- `"full"`:

  the hyperparameters are chosen again by the fit's own criteria, so
  their uncertainty enters the interval too, and the path of a kinked
  penalty is walked again, warm started from the fit's coefficients. A
  cross-validated criterion needs the data's folds and is rejected.

Every refit starts at the optimum that the fit found: its coefficients,
including those that the fit estimated on its criterion, and its
hyperparameters. The bootstrap is rejected for a model with a block that
moves with its coefficients
([`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html))
and for a response that is not one number per observation.

Each replica gives the component \\(\hat\eta^\*\_b, C^\*\_b)\\: its
predictors at the rows predicted and the covariance that the effects set
aside add, at its own hyperparameters. The interval is the mixture's,
which is the parametric bootstrap predictive distribution of Harris
(1989); the default `predictive = "averaged"` is its normal
approximation, in which the sampling distribution of the estimates is
replaced by \\\mathrm{N}(\hat \beta, \widehat{V})\\. A replica whose
refit fails is left out, and the warning reports how many were.

## References

Harris, I. R. (1989). Predictive fit for natural exponential families.
*Biometrika*, 76, 675–684.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the caller,
[`predictive_response()`](https://statmodels7.github.io/statmodels7/reference/predictive_response.md)
for the mixture.

## Examples

``` r
set.seed(2)
dd <- data.frame(x = runif(30))
dd$y <- 1 + dd$x + rnorm(30, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
predict(fit, "response", data.frame(x = 0.5), interval = "prediction",
        predictive = "bootstrap", n_boot = 20)
#>        fit        se     lower    upper
#> 1 1.547902 0.5891418 0.3767825 2.695134
```
