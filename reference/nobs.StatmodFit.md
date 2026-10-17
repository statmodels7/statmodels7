# The Number of Observations in a Fitted Model

Returns the number of observations that entered the likelihood. This
number is the sum of the prior weights, which are frequency weights.

## Usage

``` r
# S3 method for class 'StatmodFit'
nobs(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A single number, the sum of the weights: an integer where that sum is a
whole number, and a double otherwise.

## Details

A weight of two counts its row twice, so a fit with weights is counted
as the fit of the data with every row repeated as many times as its
weight. Without weights the result is the number of rows of the data,
since
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
drops no rows (a missing value in the response or in a variable of an
equation is an error).
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
the cBIC of
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)'s
`"nobs"` attribute and
[`df.residual.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/df.residual.StatmodFit.md)
read the same number.

## See also

[`df.residual.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/df.residual.StatmodFit.md),
which subtracts the effective degrees of freedom from it,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
nobs(fit)
#> [1] 40
```
