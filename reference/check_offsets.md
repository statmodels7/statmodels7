# Validate Offsets

Returns a named list with one offset per parameter, `NULL` where none
was given.

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

A named list of length `length(params)`.
