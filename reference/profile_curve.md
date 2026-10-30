# The One-Dimensional Profile of One Label

Returns the profile of one label as a function of its value alone. With
`joint = FALSE` these are the rows in which the label is held. With a
joint profile the value of the profile at each value of the label is the
maximum over the rows with that value, which is the profile over the
other labels as far as the grid resolves it. The rows without a value of
the profile are dropped, and the estimate is added with the value of the
fit.

## Usage

``` r
profile_curve(x, label, column = x@profile_column[[label]])
```

## Arguments

- x:

  A
  [StatmodProfile](https://statmodels7.github.io/statmodels7/reference/StatmodProfile.md)
  object.

- label:

  The label.

- column:

  The column of `points` that is the profile.

## Value

A data frame with the columns `value`, `profile`, `converged` and
`drop`, the likelihood ratio statistic \\2(\hat\ell - \ell_p)\\, sorted
by `value`. Its attribute `maximum` is \\\hat\ell\\.
