# Residual Degrees of Freedom of a Fitted Model

The number of observations,
[`nobs.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/nobs.StatmodFit.md),
less the effective degrees of freedom of the model.

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
freedom. Not an integer in general. Where the effective degrees of
freedom of a term could not be computed,
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
uses the number of coefficients of that term (the number of free
parameters, for a structural term) in its place, so the result is always
a number.

## Details

The count subtracted is the effective one, the sum over the terms of
their effective degrees of freedom, as
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
reports it. A penalized block has fewer effective degrees of freedom
than the number of coefficients it carries, so the result is not an
integer in general.

## See also

[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
for the per-term counts this sums,
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- sin(3 * dd$x) + rnorm(60, sd = 0.3)
fit <- statmod(y ~ s(x, bspline_smooth(k = 6)), distributions7::gaussian1_distrib(), dd)
df.residual(fit)
#> [1] 54.75431
```
