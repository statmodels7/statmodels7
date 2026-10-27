# S7 Class for the Interval of a Prediction

This class is the abstract parent of the interval objects that the
`interval` argument of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
takes:
[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md)
and
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md).
An interval object selects the uncertainty that the interval describes.

## Usage

``` r
interval_type(label = character(0))
```

## Arguments

- label:

  A short name, printed by
  [`print.interval_type()`](https://statmodels7.github.io/statmodels7/reference/print.interval_type.md).
  The constructor of each subclass sets it, since the class itself
  cannot be built.

## Value

Nothing: the class is abstract and the constructor signals an error.
Used as a class object it has one property, `label`, which every
subclass inherits.

## Details

The class is abstract and cannot be instantiated. Calling
`interval_type()` signals an error. The class makes a membership test
possible with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).
The three shipped intervals are compared in
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md).

## See also

[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md)
for the three intervals side by side,
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the consumer.

## Examples

``` r
S7::S7_inherits(confidence_interval(), interval_type)
#> [1] TRUE
S7::S7_inherits(prediction_interval(), interval_type)
#> [1] TRUE
try(interval_type("my own"))
#> Error in S7::new_object(S7::S7_object(), label = label) : 
#>   Can't construct an object from abstract class <interval_type>
```
