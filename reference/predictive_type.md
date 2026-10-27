# S7 Class for the Predictive Distribution of a New Response

This class is the abstract parent of the objects that
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
takes in its `predictive` argument:
[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
and
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md).
A predictive object selects how the estimation uncertainty enters the
distribution of a new response.

## Usage

``` r
predictive_type(label = character(0))
```

## Arguments

- label:

  A short name, printed by
  [`print.predictive_type()`](https://statmodels7.github.io/statmodels7/reference/print.interval_type.md).
  The constructor of each subclass sets it, since the class itself
  cannot be built.

## Value

Nothing: the class is abstract and the constructor signals an error.
Used as a class object it has one property, `label`, which every
subclass inherits.

## Details

The class is abstract and cannot be instantiated. Calling
`predictive_type()` signals an error. The three shipped distributions
are compared in
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md).

## See also

[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md)
for the three side by side,
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md)
for the consumer.

## Examples

``` r
S7::S7_inherits(pred_averaged(), predictive_type)
#> [1] TRUE
S7::S7_inherits(pred_bootstrap(n = 50), predictive_type)
#> [1] TRUE
try(predictive_type("my own"))
#> Error in S7::new_object(S7::S7_object(), label = label) : 
#>   Can't construct an object from abstract class <predictive_type>
```
