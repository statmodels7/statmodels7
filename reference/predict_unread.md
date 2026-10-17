# The Equations That a Prediction Does Not Read

Returns the distribution parameters whose equations
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
can leave unbuilt. These are the parameters other than `what` whose
variables `newdata` does not carry, where `what` names one parameter.

## Usage

``` r
predict_unread(object, what, newdata, interval = "confidence")
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what, newdata, interval:

  As in
  [`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## Value

A character vector of parameter names, possibly empty.

## Details

A prediction of one parameter reads only that parameter's equation, so
the covariates of the other equations need not be in `newdata`: on
`accel ~ s(times) | sigma ~ z`, `predict(fit, what = "mu", newdata)`
does not need `z`. Every equation is built where `what` is not the name
of a parameter (`"parameter"`, `"link"`, `"link:<parameter>"` or a
moment), where `interval` is not `"confidence"`, or where the model
carries a structural term, whose recursion reads the predictors of every
equation.
