# A Smoother of y on x

Fits [`stats::loess()`](https://rdrr.io/r/stats/loess.html) and predicts
it on 100 points of the range of `x`; where the fit fails, the means of
`y` within 20 classes of `x` of equal count are returned.

## Usage

``` r
diag_smooth(x, y)
```

## Arguments

- x, y:

  Numeric vectors.

## Value

A data frame with the columns `x` and `y`.
