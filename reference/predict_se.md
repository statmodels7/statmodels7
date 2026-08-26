# The Uncertainty of a Predicted Predictor

The standard error of each equation's linear predictor, and of the
parameter it gives, with an interval.

## Usage

``` r
predict_se(object, spec, design, ep, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model.

- spec:

  The specification the prediction is made under.

- design:

  Its design.

- ep:

  The predictors and parameters, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- level:

  The interval's level.

- ...:

  Passed to [`vcov()`](https://rdrr.io/r/stats/vcov.html).

## Value

A named list, one entry per distribution parameter, each a data frame
with `fit`, `se`, `lower` and `upper` on the link scale and the same
four on the parameter scale.

## The delta method, equation by equation

An equation's predictor is \\\eta\_{ip} = x\_{ip}'\beta_p\\, so its
variance is \\x\_{ip}' V\_{pp} x\_{ip}\\, with \\V\\ the variance of the
coefficients **as estimated**: the coordinates, since the design is
written in them, never the quantities
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
reports by default.

The equations do not mix. One equation's predictor reads that equation's
coefficients alone, whatever the covariance between the blocks.

## A term whose block moves needs no special case

Its block is the Jacobian \\\partial\eta/\partial\beta\\ by
construction, that being why a linear fit on it is a Gauss-Newton step,
so the row already is the derivative and the delta method is exact to
first order. That covers
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
and
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
including the parameters of a nonlinear term developed over covariates.

Measured against a numerical derivative of the predictor in the
estimated coefficients, which shares no arithmetic with the design row:
1.7e-12 on a parametric block, 8.9e-11 on a smooth with a random effect,
1.3e-10 on a `seg()` and 1.1e-11 on an `nl()` with a ridge.

## The interval

Built on the scale the equation is written on and mapped back through
the link, as every interval in this toolkit is. A scale rides a
logarithm, so its lower end cannot come out negative.

## A coefficient with no variance carries none forward

A discontinuous break-point term's block is a working linearization and
is held out of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
so every observation whose predictor reads it reports `NA` for its
standard error. That is the truth about such a fit, no gap in the
arithmetic.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
