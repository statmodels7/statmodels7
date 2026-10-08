# A Variance Matrix With Its Unreadable Coordinates Set to Zero

Sets to zero the rows and columns of a variance matrix whose diagonal is
not finite, which is how
[`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)
marks a coordinate whose curvature it could not read.

## Usage

``` r
zero_unread(V)
```

## Arguments

- V:

  A variance matrix.

## Value

`V` with those rows and columns set to zero.
