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

  Passed to [`vcov()`](https://rdrr.io/r/stats/vcov.html), which also
  receives `readable = FALSE` from this function.

## Value

A named list, one entry per distribution parameter, each a data frame
with the columns `eta`, `se_eta`, `eta_lower` and `eta_upper` on the
link scale and `fit`, `se`, `lower` and `upper` on the parameter scale.

## The delta method, equation by equation

An equation's predictor is \\\eta\_{ip} = x\_{ip}'\beta_p\\, so its
variance is \\x\_{ip}' V\_{pp} x\_{ip}\\, with \\V\\ the variance of the
coefficients as estimated: the coordinates, since the design is written
in them, not the quantities that
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
reports by default.

The equations do not mix. One equation's predictor reads that equation's
coefficients alone, whatever the covariance between the blocks. A
predictor that reads no coefficient (an equation written as `~ 0`, or an
offset alone) is a known constant, with a standard error of zero.

## A term whose block moves needs no special case

Its block is the Jacobian \\\partial\eta/\partial\beta\\ by
construction, that being why a linear fit on it is a Gauss-Newton step,
so the row already is the derivative and the delta method is exact to
first order. That covers
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
and a smoothed
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
including the parameters of a nonlinear term developed over covariates.
A sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
is fitted as an ordinary block at its settled break-point, so the
standard error of a prediction that reads it is finite and is
conditional on that position.

## The interval

Built on the scale on which the equation is written and mapped back
through the link, as every interval in this toolkit is. A scale
parameter under a log link, for example, cannot have a negative lower
end.

## A coefficient with no variance carries none forward

An entry of the variance that is not finite blanks the standard error of
the rows whose design reaches it, and of no other row. An aliased
coefficient is such a case. The position of a sharp break-point term is
held and has no variance either, but its column of the design is zero,
so it blanks no row.

A coefficient that a lasso, an elastic net, a SCAD or an MCP holds at
exactly zero is treated as fixed: the model is selected without it, its
row of the variance is `NA`, and its column of the design is left out of
the quadratic form, so it contributes nothing. Otherwise a prediction at
new data would read every such coefficient against a non-zero covariate
and return `NA` on every row.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
