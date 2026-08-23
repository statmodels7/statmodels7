# The Distribution a Model Was Fitted With

The distributions7 object, with its links.

## Usage

``` r
# S3 method for class 'StatmodFit'
family(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A distributions7 distribution.

## Details

It is the family itself rather than a description of one, so everything
the family can do is available from a fit: its density, its derivatives,
its moments and its parameters' links.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
family(fit)@params
#> [1] "mu"    "sigma"
```
