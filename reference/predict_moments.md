# The Quantities a Fit Can Predict

The table
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
resolves a moment name against: the five names it understands, each
mapped to the distributions7 generic that computes it.

## Usage

``` r
predict_moments()
```

## Value

A named list of five functions, keyed `"mean"`, `"variance"`,
`"std_dev"`, `"skewness"` and `"kurtosis"`. Each takes a distribution
and a parameter list and returns one value per observation.

## Details

Written once as a table so that the recognized names and the functions
they call cannot disagree, and so that
[`unknown_what()`](https://statmodels7.github.io/statmodels7/reference/unknown_what.md)
can list them in its message.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
the only caller.
