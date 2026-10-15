# Generics Not Available for a Fit

Three generics of stats and the generic
[`base::plot()`](https://rdrr.io/r/base/plot.html) signal an error on a
statmod fit, each naming the function to use instead.

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

Nothing. Each method signals an error that names the function to use
instead.

## Details

[`terms()`](https://rdrr.io/r/stats/terms.html) has no single answer
where a fit has one set of terms per distribution parameter, and the
formula of a fit is not a `terms` object, because the bars separating
the equations are not `stats` syntax.
[`formula.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/formula.StatmodFit.md)
gives what was written and
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
gives what it produced.

[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) would have
to return the data frame the model was fitted to, which a fit does not
keep. A fit keeps the response, the weights, the offsets and each term's
blueprint, so that a basis, a set of levels or a set of contrasts is
reapplied to new data and not learned from it again.

[`anova()`](https://rdrr.io/r/stats/anova.html) would have to compare
models by a test, and a penalized fit whose hyperparameters were chosen
from the same data has no null distribution to compare against.
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) report the fit with the sum
of the effective degrees of freedom of its terms, which is not corrected
for the smoothing parameters having been estimated;
`summary(fit, correct = TRUE)` applies that correction.

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
