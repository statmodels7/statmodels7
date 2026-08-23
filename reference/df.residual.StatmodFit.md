# What a Fit Leaves Unspent

The observation count less the effective degrees of freedom the model
spent.

## Usage

``` r
# S3 method for class 'StatmodFit'
df.residual(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A number.

## Details

The count subtracted is the EFFECTIVE one, the trace of the model's
smoother, and not the number of coefficients: a penalized block spends
less than it carries, which is the whole reason a smoothing parameter is
estimated. It is therefore not an integer, and for a model whose degrees
of freedom could not be counted it is `NA`.

## See also

[`logLik.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md),
[`statmod_edf`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- sin(3 * dd$x) + rnorm(60, sd = 0.3)
fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), dd)
df.residual(fit)
#> [1] 54.72494
```
