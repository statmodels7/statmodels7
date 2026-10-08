# What a Fit Does Not Answer

Three generics of stats and the generic
[`base::plot()`](https://rdrr.io/r/base/plot.html) signal an error on a
statmod fit, each naming what to use instead.

## Usage

``` r
# S3 method for class 'StatmodFit'
terms(x, ...)

# S3 method for class 'StatmodFit'
model.frame(formula, ...)

# S3 method for class 'StatmodFit'
anova(object, ...)

# S3 method for class 'StatmodFit'
plot(x, ...)
```

## Arguments

- ...:

  Unused.

- object, x, formula:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

Nothing. Each method signals an error naming what to ask instead.

## Details

[`terms()`](https://rdrr.io/r/stats/terms.html) would have to report one
set of terms where a fit has one per distribution parameter, and the
formula it was written with is not a `terms` object at all, the bars
separating the equations not being `stats`' syntax.
[`formula.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/formula.StatmodFit.md)
gives what was written and
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
gives what it produced.

[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) would have
to return the fitting data, which a fit does not keep. What it keeps is
each term's blueprint, so that a basis, a set of levels or a set of
contrasts is reapplied to new data instead of being learned from it
again.

[`anova()`](https://rdrr.io/r/stats/anova.html) would have to compare
models by a test, and a penalized fit whose hyperparameters were chosen
from the same data has no null distribution to compare against.
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
`AIC` and `BIC` are what this package reports, with the effective
degrees of freedom corrected for the smoothing parameters having been
estimated.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) would have to
choose one picture for a model with an equation per parameter, any
number of terms and any number of covariates. The fitted curves come
from
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
on a grid of new data, and the residuals from
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md).

## See also

[`formula.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/formula.StatmodFit.md),
[`model.matrix.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/model.matrix.StatmodFit.md),
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
