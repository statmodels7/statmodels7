# What a Fit Does Not Answer

Three generics of stats signal an error on a statmod fit, each naming
what to ask instead.

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
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

Nothing; each method signals an error.

## Details

[`terms()`](https://rdrr.io/r/stats/terms.html) would have to report one
set of terms where a fit has one per distribution parameter, and the
formula it was written with is not a `terms` object – the bars
separating the equations are not `stats`' syntax.
[`formula()`](https://rdrr.io/r/stats/formula.html) gives what was
written and
[`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
gives what it produced.

[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) would have
to return the fitting data, which a fit does not keep: what it keeps is
each term's blueprint, so that new data is reapplied rather than
relearned.

[`anova()`](https://rdrr.io/r/stats/anova.html) would have to compare
models by a test, and a penalized fit whose hyperparameters were chosen
from the same data has no null distribution to compare against.
[`logLik.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
`AIC` and `BIC` are what this package reports, with the effective
degrees of freedom corrected for the smoothing parameters having been
estimated.

## See also

[`formula.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/formula.StatmodFit.md),
[`model.matrix.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/model.matrix.StatmodFit.md),
[`logLik.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
