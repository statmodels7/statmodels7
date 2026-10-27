# Pointwise Band of the Normal Order Statistics

For a sample of size `n` from the standard normal distribution, returns
the expected plotting positions \\\Phi^{-1}(p_i)\\, with \\p_i\\ from
[`stats::ppoints()`](https://rdrr.io/r/stats/ppoints.html), and the
pointwise interval of level `level` of the order statistics. The
\\i\\-th order statistic of \\n\\ uniform variables is \\\text{Beta}(i,
n + 1 - i)\\, and the limits are the quantiles of that law mapped by
\\\Phi^{-1}\\.

## Usage

``` r
diag_order_band(n, level)
```

## Arguments

- n:

  The sample size.

- level:

  The level.

## Value

A data frame with the columns `x`, `lower` and `upper`.
