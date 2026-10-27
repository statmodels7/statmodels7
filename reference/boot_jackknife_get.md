# The Stored Refits Without One Unit

Returns the matrix of
[`bootstrap_jackknife()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_jackknife.md)
for the entry that a
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
object selects, computing it the first time and storing it in the entry.

## Usage

``` r
boot_jackknife_get(object, vt, interval)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  The resolved `BootstrapVcov` object.

- interval:

  The
  [bootstrap_interval](https://statmodels7.github.io/statmodels7/reference/bootstrap_interval.md),
  whose scheme may carry the data frame of a cluster scheme, which the
  stored entry does not.

## Value

The matrix of
[`bootstrap_jackknife()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_jackknife.md),
with the columns of the replicates of the entry.

## Details

The refits do not depend on the number of replicates of the entry, so
the stored matrix stays when the entry is extended. A message says how
many refits were made. A warning is signaled when some refits failed;
they are left out of the acceleration.
