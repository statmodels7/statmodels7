# The Quantities That a Fit Can Predict

Returns the table against which
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
resolves a moment name. The table maps each recognized name to the
distributions7 generic that computes it.

## Usage

``` r
predict_moments()
```

## Value

A named list of functions, keyed `"mean"`, `"variance"`, `"std_dev"`,
`"skewness"` and `"kurtosis"`. Each takes a distribution and a parameter
list and returns one value per observation.

## Details

The table is written once so that the recognized names and the functions
that they call cannot disagree, and so that
[`unknown_what()`](https://statmodels7.github.io/statmodels7/reference/unknown_what.md)
can list the names in its message.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the only caller.
