# The Number of Observations a Model Was Fitted To

The row count of the fitting data.

## Usage

``` r
# S3 method for class 'StatmodFit'
nobs(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

An integer.

## See also

[`logLik.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
nobs(fit)
#> [1] 40
```
