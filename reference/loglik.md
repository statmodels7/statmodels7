# The Model as a Function of Parameters and Data

Turn a fitted model back into a function. `loglik()`, `gradient()` and
`hessian()` evaluate the model's log-likelihood and its first two
derivatives at coefficients and data of the caller's choosing, so a fit
can be profiled, differenced or handed to a user-written optimizer.

If only the fit is given, these functions evaluate the model at the
fitted coefficients and on the fitting data.

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

  A named list of coefficient vectors, each as long as the design of the
  equation that it names. A distribution parameter that the list omits
  keeps its fitted coefficients, so a list may change one equation and
  leave the others at the fit. Defaults to the fitted coefficients. A
  `par` that is not a named list, a name that is not a distribution
  parameter, or a vector of the wrong length is an error. The own
  parameters of a structural term are not among the arguments.

- data:

  A data frame with the columns that the model names. Defaults to the
  data used to fit the model. New data must carry the response, since a
  likelihood needs one.

- ...:

  Passed to methods. The shipped methods do not read it.

- expected:

  `hessian()` only: `TRUE` for the expected information, `FALSE` (the
  default) for the observed Hessian.

## Value

`loglik()` gives a single number. `gradient()` gives a named list of
numeric vectors, one per distribution parameter. `hessian()` gives a
symmetric `p x p` matrix over the stacked regression coefficients, where
`p` is their total count. The own parameters of a structural term are
not among them and have no entry in the gradient or in the Hessian.

## The name

The generic is `loglik`, not `logLik`. R's
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) returns the
log-likelihood of a fitted model and carries `df` and `nobs` with it,
and overloading it would give one name two behaviors. Both functions are
available, and `loglik(fit)` with no further arguments equals
`logLik(fit)`, which checks that the callable route and the fitting
route evaluate the same model.

## Generics instead of closures

A closure captures its environment, which here means the data. A fit
carrying one would hold the frame twice and would keep a stale copy
after the data changed.

The generics rebuild from the specification that the fit keeps whole,
running each term's blueprint against `data` by the same path as
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).
A factor's levels, a spline's knots and a basis reparametrization are
therefore reapplied to the new rows and not estimated again.

## Likelihood without the penalty

All three functions evaluate the likelihood without the penalty. The
penalized objective at given coefficients is a different quantity, and
the hyperparameters are not arguments of these functions.

## See also

[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
for the log-likelihood of the fit with its degrees of freedom,
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

# At the fitted coefficients it equals the conditional log-likelihood.
all.equal(loglik(fit), as.numeric(logLik(fit)))
#> [1] TRUE

# Far from the fit it is lower, and the gradient gives the direction of
# increase.
loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#> [1] -88.88139
gradient(fit, par = list(mu = c(0, 0), sigma = 0))
#> $mu
#> [1] 61.90424 34.21924
#> 
#> $sigma
#> [1] 64.2477
#> 

# The gradient with respect to the mu coefficients vanishes at the fit.
# The sigma coordinate is estimated on the REML criterion, so its
# component is minus the number of mu coefficients (here -2).
gradient(fit)
#> $mu
#> [1] 1.421085e-14 7.105427e-15
#> 
#> $sigma
#> [1] -2
#> 

# The observed information is positive definite at the fit.
eigen(-hessian(fit), only.values = TRUE)$values
#> [1] 294.47910  76.00000  13.73987
```
