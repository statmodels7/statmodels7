# Print a Bootstrap Interval

Prints one line, the label of the interval and the settings of the
bootstrap that it reads.

## Arguments

- x:

  A
  [bootstrap_interval](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(percentile_interval(n = 500))
#> <interval> percentile bootstrap interval (500 replicates, resample as stored, refit as stored)
```
