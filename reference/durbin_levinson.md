# Partial Autocorrelations by the Durbin-Levinson Recursion

Turns the autocorrelations \\r_1, \dots, r_K\\ into the partial
autocorrelations \\\phi\_{11}, \dots, \phi\_{KK}\\.

## Usage

``` r
durbin_levinson(r)
```

## Arguments

- r:

  A numeric vector of autocorrelations at the lags 1 to K.

## Value

A numeric vector of the same length.
