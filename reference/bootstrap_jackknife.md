# The Refits Without One Observation or One Cluster

Refits the model once for each unit with the weights of its rows set to
zero, and returns the coefficients of the refits.

## Usage

``` r
bootstrap_jackknife(object, resample, refit)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  The
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  of the entry.

- refit:

  `"coefficients"` or `"full"`.

## Value

A matrix with a row for each unit and a column for each coordinate of
the replicates of the entry, with the attribute `failed`.

## Details

The units are the observations with a positive weight, or the clusters
of
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
that contain one. The refit is the one of the stored entry (`refit`),
and the replicate is checked as the engine checks it. The random number
stream of the session is not moved. A refit that fails leaves a row of
`NA`, and the attribute `failed` holds the number of such rows.
