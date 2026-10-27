# The Fit as the Workers Receive It

Returns the fit without the optimizer of its outer search.

## Usage

``` r
boot_light(object)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

The
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
with `methods$search` set to `NULL`.

## Details

The closures of the engine hold the fit, and
[`worker_map()`](https://statmodels7.github.io/statmodels7/reference/worker_map.md)
serializes them for every worker. The optimizer that the outer search
stores in `object@methods$search` keeps the environment of that search,
which holds the design and the inner fit, and it is not read by a refit.
Dropping it leaves the properties that
[`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md),
[`boot_read_mode()`](https://statmodels7.github.io/statmodels7/reference/boot_read_mode.md)
and the draws read.
