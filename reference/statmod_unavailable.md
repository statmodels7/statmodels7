# Generics Not Available for a Fit

Three generics of stats signal an error on a statmod fit, and each error
names the function to use instead.

## Usage

``` r
# S3 method for class 'StatmodFit'
terms(x, ...)

# S3 method for class 'StatmodFit'
model.frame(formula, ...)

# S3 method for class 'StatmodFit'
anova(object, ...)
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
to return the data frame used to fit the model, which a fit does not
keep. A fit keeps the response, the weights, the offsets and each term's
blueprint, so that a basis, a set of levels or a set of contrasts is
reapplied to new data and not learned from it again.

[`anova()`](https://rdrr.io/r/stats/anova.html) would have to compare
models by a test, and a penalized fit whose hyperparameters were chosen
from the same data has no null distribution against which to compare.
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) report the fit with the sum
of the effective degrees of freedom of its terms, which is not corrected
for the estimation of the smoothing parameters;
`summary(fit, correct = TRUE)` applies that correction.

## See also

[`formula.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/formula.StatmodFit.md),
[`model.matrix.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/model.matrix.StatmodFit.md),
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
