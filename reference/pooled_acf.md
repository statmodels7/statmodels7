# Pooled Autocorrelation of Several Series

Computes \\r_1, \dots, r_K\\ of
[`plot_acf()`](https://statmodels7.github.io/statmodels7/reference/plot_acf.md),
with the products of the pairs within each series in the numerator and
the sum of squares of all the observations in the denominator, both
about the overall mean. A pair with a missing value is left out of the
numerator.

## Usage

``` r
pooled_acf(e, series, lag_max)
```

## Arguments

- e:

  A numeric vector.

- series:

  A list of integer vectors, the positions of each series in time order.

- lag_max:

  The largest lag.

## Value

A numeric vector of length `lag_max`.
