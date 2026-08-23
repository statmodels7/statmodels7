# Print a Model Summary

The call, then each distribution parameter's blocks, then the degrees of
freedom, the criteria and the notes.

## Usage

``` r
# S3 method for class 'StatmodSummary'
print(x, digits = 4L, notes = FALSE, ...)
```

## Arguments

- x:

  A
  [`StatmodSummary`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md).

- digits:

  Significant digits in the tables.

- notes:

  Whether to print the qualifications the numbers carry. `FALSE` by
  default, when the foot says how many there are: they state conventions
  rather than facts of the fit, so they read the same under every model.
  They are on the summary's `notes` property either way.

- ...:

  Unused.

## Value

`x`, invisibly.

## See also

[`summary.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
