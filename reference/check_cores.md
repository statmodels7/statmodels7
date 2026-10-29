# Check the Number of Processes

Returns `cores` as an integer if it is a single whole number of at least
1, and signals an error otherwise. `cores` is the number of separate R
processes over which a function distributes independent fits; it is not
the `threads` of a fit, which counts threads that share memory inside
one process.

## Usage

``` r
check_cores(cores)
```

## Arguments

- cores:

  The value of a `cores` argument.

## Value

A single integer, at least 1.

## Examples

``` r
statmodels7:::check_cores(2)
#> [1] 2
try(statmodels7:::check_cores(0))
#> Error : 'cores' must be a single whole number, at least 1.
```
