# Check the Interval Argument of predict()

Validates the `interval` argument of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
and returns the kind of interval and, for a prediction interval, the
kind and the settings of the predictive distribution, as the strings
that the body of the prediction reads.

## Usage

``` r
predict_interval_spec(interval)
```

## Arguments

- interval:

  The argument as passed.

## Value

A list with `interval` (`"confidence"`, `"group"` or `"prediction"`),
`predictive` (`"averaged"`, `"plugin"` or `"bootstrap"`), `n_boot` and
`boot_refit`.
