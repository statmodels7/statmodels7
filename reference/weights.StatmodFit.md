# The Prior Weights of a Fitted Model

Returns the weight with which each observation entered the likelihood.

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

A numeric vector with one entry per row of the data, the prior weights
as supplied. The vector is all ones if no weights were given. The
weights enter as given and are never normalized, so a weight of two
counts an observation twice.

## Details

The weights multiply the log-likelihood of each observation, so they act
as frequency weights: a constant weight \\c\\ on every row is \\c\\
copies of the data. The restricted likelihood then counts \\n c\\
observations against the \\p\\ coefficients, so the REML estimate of a
dispersion changes under a constant rescaling of the weights, while the
ML estimate does not.

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

# The weights are frequencies, so nobs() is their sum.
sum(weights(fit)) == nobs(fit)
#> [1] TRUE
```
