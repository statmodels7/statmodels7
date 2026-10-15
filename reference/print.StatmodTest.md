# Print a Test

The statistic, its degrees of freedom and the p-value on one line, under
a title naming the test, in the layout R uses for a hypothesis test. A
note follows when the restricted refit stopped above its mode, by more
than
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
or when its mode error could not be read.

## Usage

``` r
# S3 method for class 'StatmodTest'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  A
  [`StatmodTest()`](https://statmodels7.github.io/statmodels7/reference/StatmodTest-class.md).

- digits:

  Significant digits.

- ...:

  Unused.

## Value

`x`, invisibly.
