# The Variance of Each Equation's Working Response

One number per distribution parameter: the weighted variance of \\z_i =
\eta_p + U\_{ip} / I\_{ip} - o\_{ip}\\, the working response at the
intercept-only fit less the offset, as
[`statmod_hyper_scale_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_scale_start.md)
defines it.

## Usage

``` r
working_response_variance(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A named list, one entry per parameter; an entry is `NULL` where the
intercept-only fit, the score or the information is not available.
