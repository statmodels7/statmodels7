# Print a Resampling Scheme

Prints one line, the name of the scheme and, for a cluster scheme, its
formula.

## Arguments

- x:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(cluster_resample(~ id))
#> <resample> cluster (~id)
```
