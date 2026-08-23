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
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A numeric vector.

## Details

A VECTOR rather than a number, because the whole point of the framework
is that a scale may be modelled: a single residual standard deviation
exists only where the scale's equation is an intercept, and returning
its first value would silently answer a different question everywhere
else.

What is returned is the standard deviation of the response under the
fitted distribution, through
[`std_dev`](https://statmodels7.github.io/distributions7/reference/std_dev.html),
and not whichever parameter happens to be spelled `sigma`: for a Gamma
written by its mean and dispersion the two are different quantities. A
family with no second moment signals an error rather than reporting one.

## See also

[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`fitted.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
head(sigma(fit))
#> [1] 0.2220141 0.2303400 0.2468724 0.2771821 0.2171743 0.2762440
```
