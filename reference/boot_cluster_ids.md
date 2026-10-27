# The Groups of a Cluster Bootstrap

Evaluates the variable of
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
on the data of the fit.

## Usage

``` r
boot_cluster_ids(object, resample, spec0)

boot_cluster_data(object, resample, spec0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  A
  [ClusterResample](https://statmodels7.github.io/statmodels7/reference/ParametricResample.md).

- spec0:

  Its specification.

## Value

A vector with one value per observation.

## Details

The data frame is the `data` of the scheme if given, and otherwise the
`data` of the call that built the fit, evaluated in the environment of
the model formula. The fit does not keep its data, so a `data` argument
that no longer evaluates is an error that says to pass `data`.
