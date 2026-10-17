# The Shape in Which a Prediction With Its Uncertainty Is Returned

Reduces the rows that
[`predict_se()`](https://statmodels7.github.io/statmodels7/reference/predict_se.md)
computed to the quantity requested.

## Usage

``` r
se_answer(su, what, params, spec)
```

## Arguments

- su:

  The per-parameter tables.

- what:

  The quantity requested, as in
  [`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

- params:

  The distribution's parameters.

- spec:

  The specification, for the message.

## Value

A data frame, or a named list of them.

## Details

The accepted values of `what` are those of a prediction without
uncertainty, minus the moments: a moment's delta method needs its
derivative in every parameter, which distributions7 does not provide as
a generic, so a moment signals an error.
