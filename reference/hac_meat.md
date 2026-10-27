# The HAC Meat

Computes \\\Gamma_0 + \sum\_{l=1}^{L} k(l/(L+1))(\Gamma_l + \Gamma_l')\\
from the scores and the series.

## Usage

``` r
hac_meat(U, series, kernel = "bartlett", bandwidth = NULL)
```

## Arguments

- U:

  The matrix of scores.

- series:

  The rows of each independent series, in time order.

- kernel:

  The name of the kernel.

- bandwidth:

  The largest lag, or `NULL` for the default.

## Value

A list with the symmetric `meat` and the `bandwidth` used.
