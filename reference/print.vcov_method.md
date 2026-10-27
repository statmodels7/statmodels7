# Print a Variance Convention

Prints one line, the name of the convention and, for the three
conventions that read the curvature, the information that they invert.

## Arguments

- x:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(unconditional_vcov())
#> <vcov> unconditional
print(bayesian_vcov(information = distributions7::expected_info()))
#> <vcov> bayesian (expected information)
```
