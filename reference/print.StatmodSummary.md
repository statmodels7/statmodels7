# Print a Model Summary

The call, then the covariance blocks shared between equations, then each
distribution parameter's blocks, then the degrees of freedom, the
criteria and the notes.

Below the criteria the point is described by dimensionless numbers from
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md):
for the inner fit, and for the outer criterion where there is one, the
largest gradient in conditional standard errors and the smallest
eigenvalue of the equilibrated curvature
([`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)),
and for the coefficients that a kinked penalty has set to zero the
largest ratio of the pull of the likelihood to the kink
([`zero_readings()`](https://statmodels7.github.io/statmodels7/reference/zero_readings.md)).
No verdict is printed. How each is read is on the page of
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
under "Reading the foot". The line of the observations also gives the
sum of the prior weights where they are not all 1.

## Usage

``` r
# S3 method for class 'StatmodSummary'
print(x, digits = 4L, notes = FALSE, max_coef = NULL, ...)
```

## Arguments

- x:

  A
  [`StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md).

- digits:

  Significant digits in the tables.

- notes:

  Whether to print the notes that qualify the numbers. With the default
  `FALSE` the foot gives only the number of notes. Most of them state
  conventions of the model, and some name the coefficients or the values
  that the fit reached, so they can differ between models. The notes are
  on the summary's `notes` property in either case.

- max_coef:

  How many coefficient rows a block shows, `Inf` or `NA` for all of
  them. `NULL`, the default, reads what
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  was told; failing that the option `statmodels7.summary_max_coef`, and
  failing that 10. A hyperparameter row is shown whatever `max_coef` is:
  it governs the coefficients under it, and every one of them is
  conditional on the value that it reached. A block short enough to fit
  in twelve rows is never abridged, so raising `max_coef` changes
  nothing for a parametric block of ordinary size.

- ...:

  Unused.

## Value

`x`, invisibly.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
