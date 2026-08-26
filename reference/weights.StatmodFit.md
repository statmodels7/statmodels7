# The Prior Weights of a Fitted Model

The weights each observation entered the likelihood with.

## Usage

``` r
# S3 method for class 'StatmodFit'
weights(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A numeric vector of length `nobs(object)`, the prior weights as
supplied. All ones when none were given. They enter as given and are
never normalized, so a weight of two counts an observation twice.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`nobs.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/nobs.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
head(weights(fit))
#> [1] 1 1 1 1 1 1

# An unweighted fit carries ones, so the weights sum to the row count.
sum(weights(fit)) == nobs(fit)
#> [1] TRUE
```
