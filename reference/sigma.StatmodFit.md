# The Standard Deviation Implied by a Fitted Model

Returns the fitted standard deviation at each observation, where the
family has one.

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
response under the fitted distribution at each observation. It is
constant when no equation models a scale and the variance of the family
does not depend on the location; for a Poisson family, for example, it
varies with the fitted mean.

## Details

The result is a vector with one entry per observation, because a scale
may be modeled. A single residual standard deviation exists only where
the equation of the scale is an intercept, and the first value of the
vector would be a different quantity in any other model.

The value is the standard deviation of the response under the fitted
distribution, computed through
[`distributions7::std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.html)
from the fitted parameters. It is not the parameter that happens to be
named `sigma`: for a Gamma written by its mean and dispersion, the two
are different quantities. Where the second moment does not exist at the
fitted parameters (the Cauchy distribution, for example), the entries
are `NaN`.

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
#> [1] 0.2282563 0.2366059 0.2531638 0.2834521 0.2233994 0.2825159
```
