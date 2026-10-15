# The Error for a Cross-Validation With Nothing to Select

Signals the error that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
gives when
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) is
passed as `outer_criterion` to a model in which no penalty with a kink
has a hyperparameter left to estimate.

## Usage

``` r
stop_cv_nothing()
```

## Value

Does not return.
