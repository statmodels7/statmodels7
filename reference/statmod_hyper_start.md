# The Hyperparameters a Specification Starts From

One entry per penalized term, at the midpoint of its bounds – the probe
rule modelterms7 already uses when it reads a penalty's kinks.

## Usage

``` r
statmod_hyper_start(spec, design = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A named list, one entry per parameter.
