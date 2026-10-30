# The Points That a Path Scored

The Points That a Path Scored

## Usage

``` r
profile_path_record(fit, r, link)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- r:

  The row of the hyperparameter in the table of
  [`profile_hyper_table()`](https://statmodels7.github.io/statmodels7/reference/profile_hyper_table.md).

- link:

  The link of the hyperparameter.

## Value

A data frame with the columns `value`, on the scale of the link,
`criterion` and `se`, sorted by `value`; it has no rows where the fit
kept no record of the path.
