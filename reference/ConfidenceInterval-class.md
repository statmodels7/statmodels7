# S7 Classes for the Shipped Intervals and Predictive Distributions

These are the classes that
[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[`pred_averaged()`](https://statmodels7.github.io/statmodels7/reference/pred_averaged.md),
[`pred_plugin()`](https://statmodels7.github.io/statmodels7/reference/pred_plugin.md)
and
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
instantiate. `ConfidenceInterval`, `GroupInterval` and
`PredictionInterval` inherit from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md);
`PredAveraged`, `PredPlugin` and `PredBootstrap` inherit from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md).
`PredictionInterval` adds the property `predictive`, and `PredBootstrap`
adds `n` and `refit`. Use the constructors, which validate; these raw
ones do not.

## Usage

``` r
ConfidenceInterval(label = character(0))

GroupInterval(label = character(0))

PredictionInterval(label = character(0), predictive = predictive_type())

PredAveraged(label = character(0))

PredPlugin(label = character(0))

PredBootstrap(label = character(0), n = integer(0), refit = character(0))
```

## Arguments

- label:

  A single string, the name that the object prints.

- predictive:

  A
  [`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md).
  Only `PredictionInterval` carries this property.

- n, refit:

  The number of replicas and the refit mode. Only `PredBootstrap`
  carries these properties.

## Value

An S7 object inheriting from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md)
or from
[`predictive_type()`](https://statmodels7.github.io/statmodels7/reference/predictive_type.md),
with the properties listed above.

## See also

[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md),
[predictive_distributions](https://statmodels7.github.io/statmodels7/reference/predictive_distributions.md)
for the constructors side by side.
