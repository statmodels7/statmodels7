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
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A single number, `nobs(object)` less the total effective degrees of
freedom. Not an integer in general, and `NA` for a model whose degrees
of freedom could not be counted.

## Details

The count subtracted is the effective one, the trace of the model's
smoother. A penalized block spends less than the number of coefficients
it carries, which is the whole reason a smoothing parameter is estimated
at all, so the result is not an integer.

## See also

[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
for the per-term counts this sums,
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- sin(3 * dd$x) + rnorm(60, sd = 0.3)
fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), dd)
df.residual(fit)
#> [1] 54.72494
```
