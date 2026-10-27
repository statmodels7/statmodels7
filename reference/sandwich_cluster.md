# The Cluster Index of a Sandwich Variance

Reads the cluster variable of a
[SandwichMethod](https://statmodels7.github.io/statmodels7/reference/SandwichMethod.md)
object from its data or from the data of the call that made the fit,
checks it, and checks that every group of a random term lies in one
cluster.

## Usage

``` r
sandwich_cluster(fit, design, sw)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- sw:

  The
  [SandwichMethod](https://statmodels7.github.io/statmodels7/reference/SandwichMethod.md)
  object.

## Value

An integer vector of cluster numbers, one for each observation.
