# The Data Frame Against Which a Simulation Runs

Resolves `data` and `n` into one data frame of covariates.

## Usage

``` r
rstatmod_data(data, n, covariates = NULL)
```

## Arguments

- data:

  A data frame or `NULL`.

- n:

  A row count or `NULL`.

- covariates:

  The generators, or `NULL`. Where there are any, their columns are
  written in.

## Value

A data frame.

## Details

A model may have no covariates at all (a pure time series is an
example), and then a row count is all that a simulation needs, so `data`
and `n` are alternatives and one of them must be given. Where both are
given they must agree. The function checks this and does not resolve a
disagreement in favor of one of them, since a disagreement points to a
mistake of the caller. `n` must be a positive whole number.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
