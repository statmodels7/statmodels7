# Print a Model Summary

The call, then each distribution parameter's blocks, then the degrees of
freedom, the criteria and the notes.

## Usage

``` r
# S3 method for class 'StatmodSummary'
print(x, digits = 4L, notes = FALSE, n = NULL, ...)
```

## Arguments

- x:

  A
  [`StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md).

- digits:

  Significant digits in the tables.

- notes:

  Whether to print the qualifications the numbers carry. `FALSE` by
  default, when the foot says how many there are: they state
  conventions, never facts of the fit, so they read the same under every
  model. They are on the summary's `notes` property either way.

- n:

  How many coefficient rows a block shows, `Inf` or `NA` for all of
  them. Defaults to the option `statmodels7.summary_rows`, and to 10
  where that is unset. A hyperparameter row is shown whatever `n` is: it
  governs the coefficients under it, and every one of them is
  conditional on the value it reached. A block short enough to fit in
  twelve rows is never abridged, so raising `n` changes nothing for a
  parametric block of ordinary size.

- ...:

  Unused.

## Value

`x`, invisibly.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
