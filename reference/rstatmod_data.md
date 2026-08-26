# The Data Frame a Simulation Runs Against

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

A model may have no covariates at all – a pure time series is the case –
and then a row count is the whole of what a simulation needs, so the two
arguments are alternatives rather than one being compulsory. Where both
are given they must agree, which is checked rather than resolved by
preferring one: a caller who wrote both and got them wrong wants to
know.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
