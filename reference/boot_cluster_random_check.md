# Reject a Random-Effect Term That the Cluster Bootstrap Cannot Regroup

Signals an error for a fit whose
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
terms the cluster bootstrap of
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
cannot rebuild on a resampled data set.

## Usage

``` r
boot_cluster_random_check(spec0)
```

## Arguments

- spec0:

  The specification of the fit.

## Value

`NULL`, invisibly.

## Details

A replicate rebuilds every random-effect term on its own groups, with
the `build_info` of the fit for the within-group design. A term that
shares a covariance through a label, a smooth random effect
([`modelterms7::RandomSmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomSmoothTerm.html))
and a random effect inside the subformula of another term are not
rebuilt in this version, and the fit is rejected.
