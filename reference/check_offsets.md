# Validate Offsets

Returns a named list with one offset per parameter, `NULL` for a
parameter whose offset was not given.

## Usage

``` r
check_offsets(offsets, params, n)
```

## Arguments

- offsets:

  A named list, or `NULL`.

- params:

  The parameter names.

- n:

  The number of observations.

## Value

A named list of length `length(params)`, and each entry is a numeric
vector of length `n` or `NULL`.

## Details

An entry must be numeric, of length 1 or `n`, and finite everywhere, as
an offset written in the formula must be; the names must be parameters
of the distribution.
