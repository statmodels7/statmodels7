# The Model as a Function of Parameters and Data

Turn a fitted model back into a function. `loglik()`, `gradient()` and
`hessian()` evaluate the model's log-likelihood and its first two
derivatives at coefficients and data of the caller's choosing, so a fit
can be profiled, differenced or handed to an optimizer of one's own.

With no arguments beyond the fit they read at the fitted coefficients
and the fitting data.

## Usage

``` r
loglik(object, par = NULL, data = NULL, ...)

gradient(object, par = NULL, data = NULL, ...)

hessian(object, par = NULL, data = NULL, expected = FALSE, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- par:

  A named list of coefficient vectors, one per distribution parameter,
  each as long as that equation's design is wide. Defaults to the fitted
  coefficients. Validated against the design, so a wrong length is an
  error, never a recycling.

- data:

  A data frame with the columns the model names. Defaults to the data
  the model was fitted to. New data must carry the response, since a
  likelihood needs one.

- ...:

  Passed to methods. No shipped method reads it.

- expected:

  Whether the expected information is wanted.

## Value

`loglik()` gives a single number. `gradient()` gives a named list of
numeric vectors, one per distribution parameter. `hessian()` gives a
symmetric `p x p` matrix over the stacked coefficients, `p` being their
total count.

## The name

`loglik`, not `logLik`. R's
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) returns the
maximized value of the fitted model and carries `df` and `nobs` with it,
and overloading it would give one name two behaviors. Both exist here,
and `loglik(fit)` with no further arguments equals `logLik(fit)` to the
last digit, which is the cheapest check that the callable route and the
fitting route are the same model.

## Generics, not closures on the fit

A closure captures its environment, which here means the data. A fit
carrying one would hold the frame twice and would keep a stale copy
after the data changed.

These rebuild from the specification the fit keeps whole, running each
term's blueprint against `data` by the same path
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
takes. A factor's levels, a spline's knots and a basis reparametrization
are therefore reapplied, never relearned.

## No penalty enters

All three are the likelihood alone. A model's penalized objective at
given coefficients is not what a caller profiling a likelihood wants,
and the hyperparameters are not arguments here.

## See also

[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
for the maximized value with its degrees of freedom,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the variance matrix,
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the same reapplication on new data.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)

# At the fitted coefficients it is the maximized value.
all.equal(loglik(fit), as.numeric(logLik(fit)))
#> [1] TRUE

# Anywhere else it is lower, and the gradient says which way to go.
loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#> [1] -88.88139
gradient(fit, par = list(mu = c(0, 0), sigma = 0))
#> $mu
#> [1] 61.90424 34.21924
#> 
#> $sigma
#> [1] 64.2477
#> 

# At the optimum the gradient vanishes.
max(abs(unlist(gradient(fit))))
#> [1] 2.815559e-06

# And the observed information is positive definite there.
eigen(-hessian(fit), only.values = TRUE)$values
#> [1] 309.97798  79.99999  14.46302
```
