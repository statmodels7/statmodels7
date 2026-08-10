# The Model as a Function of Parameters and Data

`loglik()`, `gradient()` and `hessian()` evaluate a fitted model at
parameters and data of the caller's choosing.

## Usage

``` r
loglik(object, par = NULL, data = NULL, ...)

gradient(object, par = NULL, data = NULL, ...)

hessian(object, par = NULL, data = NULL, expected = FALSE, ...)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- par:

  A named list of coefficient vectors, one per distribution parameter.
  Defaults to the fitted ones.

- data:

  A data frame. Defaults to the data the model was fitted to.

- ...:

  Passed to methods.

- expected:

  Whether the expected information is wanted.

## Value

A number for `loglik`, a named list of vectors for `gradient`, a matrix
for `hessian`.

## Details

The name is `loglik` and not `logLik`: R's
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the maximized
value of the fitted model and carries `df` and `nobs`, and overloading
it would give one name two behaviours. Both exist, and `loglik(fit)`
with no arguments must equal `logLik(fit)` to the last digit, which is
the cheapest check that the callable route and the fitting route are the
same model.

They are generics rather than closures stored in the fit. A closure
captures its environment, which means the data: a fit would then carry
the frame twice and keep a stale copy after the data changed. These
rebuild from the specification the fit keeps whole, running the terms'
blueprints against new data by the same path
[`predict()`](https://rdrr.io/r/stats/predict.html) takes.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
loglik(fit)
#> [1] -20.76686
loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#> [1] -88.88139
```
