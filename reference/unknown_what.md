# The Message for an Unrecognized Prediction Target

Builds the error
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
signals when `what` names nothing it can compute, listing this family's
own parameter names, the five moments and the two collective targets.

## Usage

``` r
unknown_what(what, params)
```

## Arguments

- what:

  What was asked for, for the message.

- params:

  The family's parameter names, in the family's order.

## Value

A single string, ready for [`stop()`](https://rdrr.io/r/base/stop.html).

## Details

The family's parameters are listed by name rather than described, since
they differ from family to family and are the commonest thing a caller
means. A data frame passed where `what` belongs, which is the mistake
[`stats::predict()`](https://rdrr.io/r/stats/predict.html)'s argument
order invites, is recognized and named separately.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the caller,
[`predict_moments()`](https://statmodels7.github.io/statmodels7/reference/predict_moments.md)
for the moment names listed.
