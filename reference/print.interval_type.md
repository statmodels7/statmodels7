# Print an Interval or a Predictive Distribution

Prints a one-line description of an
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md)
or a
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md)
object: its label, and for
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
the predictive distribution that it holds.

## Arguments

- x:

  An
  [`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md)
  or
  [`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(confidence_interval())
#> <interval> confidence
print(prediction_interval(predictive = pred_bootstrap(n = 50)))
#> <interval> prediction (predictive: bootstrap, 50 replicas, refit coefficients)
```
