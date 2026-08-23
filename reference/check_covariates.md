# The Covariate Generators of a Simulation

Validates `covariates` and returns it, or `NULL`.

## Usage

``` r
check_covariates(covariates)
```

## Arguments

- covariates:

  A named list, or `NULL`.

## Value

The list, or `NULL`.

## Details

Each entry is a function of the observation count, as an entry of `par`
may be, so the two arguments read alike. A value that is not a function
is rejected rather than recycled: a constant column is what `data` is
for, and the whole point of this argument is that the covariates are
drawn afresh at every replicate.

## See also

[`rstatmod`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
