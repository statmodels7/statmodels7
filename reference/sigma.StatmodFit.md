# The Standard Deviation a Fitted Model Implies

The fitted standard deviation at each observation, where the family has
one.

## Usage

``` r
# S3 method for class 'StatmodFit'
sigma(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A numeric vector of length `nobs(object)`, the standard deviation of the
response under the fitted distribution at each observation. Constant
when no equation models a scale.

## Details

The result is a vector, one entry per observation, and that is the whole
point of the framework: a scale may be modeled. A single residual
standard deviation exists only where the scale's equation is an
intercept, and returning its first value would answer a different
question everywhere else without saying so.

What comes back is the standard deviation of the response under the
fitted distribution, through
[`distributions7::std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.html).
It is not whichever parameter happens to be spelled `sigma`: for a Gamma
written by its mean and dispersion those are two different quantities. A
family with no second moment signals an error.

## See also

[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
for the fitted parameters themselves,
[`distributions7::std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.html)
for the quantity

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
head(sigma(fit))
#> [1] 0.2220141 0.2303400 0.2468724 0.2771821 0.2171743 0.2762440
```
