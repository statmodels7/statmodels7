# Resolve a Bootstrap Variance Against the Stored Entries

Returns the
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
object with `n`, `resample` and `refit` all set, after making sure that
the fit holds an entry with at least `n` replicates for those settings.

## Usage

``` r
boot_resolve(fit, vt)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  A `BootstrapVcov` object.

## Value

A `BootstrapVcov` object with no `NULL` property.

## Details

A `NULL` argument takes the value of the single entry that matches the
arguments given, and the default (`n = 200`,
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md),
`"coefficients"`) when none matches. Several matching entries are an
error. A missing entry, or one with fewer replicates than `n`, is
computed with a message, which is the only place where the engine runs
on behalf of a variance. An object that is already resolved finds its
entry and computes nothing, so
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
resolve their `vcov_type` once and pass the result on.
