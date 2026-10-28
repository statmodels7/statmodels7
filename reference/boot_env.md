# The Storage of the Bootstrap Entries of a Fit

`boot_env()` returns the environment in which a fit keeps its entries,
or `NULL` for a fit that has none (an object saved before the property
existed). `boot_fit_estimates()` returns the numbers that identify the
fit an entry belongs to, and `boot_entries()` the entries that belong to
the fit, in alphabetical order of their keys.

## Usage

``` r
boot_env(fit)

boot_fit_estimates(fit)

boot_entries(fit)

boot_match(ents, resample, refit)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ents:

  A list of entries.

- resample, refit:

  The settings to match, `NULL` for any.

## Value

The environment, a numeric vector, or a list of
[StatmodBootstrap](https://statmodels7.github.io/statmodels7/reference/StatmodBootstrap.md)
objects.
