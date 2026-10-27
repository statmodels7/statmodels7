# Check the Arguments of the Bootstrap

`check_boot_n()` returns `n` as an integer if it is one positive whole
number, and `check_boot_refit()` returns `refit` if it is
`"coefficients"` or `"full"`. Both signal an error otherwise.

## Usage

``` r
check_boot_n(n)

check_boot_refit(refit)
```

## Arguments

- n:

  The number of replicates.

- refit:

  The refit.

## Value

The checked value.
