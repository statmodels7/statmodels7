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
  [`statmod_eta`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- level:

  The interval's level.

- ...:

  Passed to [`vcov`](https://rdrr.io/r/stats/vcov.html).

## Value

A named list, one entry per distribution parameter, each a data frame
with `fit`, `se`, `lower` and `upper` on the link scale and the same
four on the parameter scale.

## Details

An equation's predictor is \\\eta\_{ip} = x\_{ip}'\beta_p\\, so its
variance is \\x\_{ip}' V\_{pp} x\_{ip}\\ with \\V\\ the variance of the
coefficients as estimated – the coordinates and not the quantities, the
design being written in them. The equations do not mix here: the
predictor of one reads that one's coefficients alone, whatever the
covariance between the blocks.

A TERM WHOSE BLOCK MOVES with its coefficients needs no special case.
Its block is the Jacobian \\\partial\eta/\partial\beta\\ by construction
– that is what makes a linear fit on it a Gauss-Newton step – so the row
is the derivative and the delta method is exact to first order. That
covers `seg`, `jseg` and `nl`, including the parameters of a nonlinear
term developed over covariates.

The interval is built on the scale the equation is written on and mapped
back through the link, as every interval in the toolkit is: a scale
rides a logarithm and its lower end cannot come out negative.

A coefficient with no variance carries none forward. A discontinuous
break-point term's block is a working linearization, held out of
[`vcov`](https://rdrr.io/r/stats/vcov.html), so every observation whose
predictor reads it has no standard error either – which is the truth
about it and not a gap in the arithmetic.

## See also

[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
