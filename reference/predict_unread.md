# The Equations a Prediction Does Not Read

The distribution parameters whose equations
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
can leave unbuilt: those other than `what` whose variables `newdata`
does not carry, where `what` names one parameter.

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

A prediction of one parameter reads that parameter's equation alone, so
the covariates of the other equations need not be in `newdata`. They
used to be required all the same, the design being rebuilt for every
equation: `predict(fit, what = "mu", newdata)` on
`accel ~ s(times) | sigma ~ z` stopped with "object 'z' not found".
Nothing is left out where every equation can be built, where `what` asks
for every parameter or a moment, or where the model carries a structural
term, whose recursion reads the predictors of every equation.
