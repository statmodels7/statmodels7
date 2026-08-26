# Confidence Intervals for a Fit

Wald intervals for the coefficients of every distribution parameter.

## Usage

``` r
# S3 method for class 'StatmodFit'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("bayesian", "frequentist"),
  readable = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- parm:

  Which coefficients: a distribution parameter's name, a vector of
  `parameter:coefficient` labels, or `NULL` for all of them.

- level:

  The confidence level.

- type:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A data frame with the parameter, the term, the coefficient, the
estimate, its standard error and the two limits.

## Details

The interval is symmetric about the estimate and needs no mapping back.
A coefficient of a linear predictor is unbounded whatever the
distribution parameter it belongs to, the link having already carried
that parameter onto the whole line, so the scale the interval is built
on is the scale the quantity lives on. What the interval does not do is
respect a bound on the parameter itself; for that, map an interval for
the predictor through the inverse link at the covariate values of
interest.

The variance comes from
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
so the same two conventions apply, and a coefficient a kinked penalty
set to zero has `NA` in place of an interval.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
confint(fit)
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.043291 0.08659784  0.8735627
#> mu:x                     mu linpar           x  2.007856 0.14678520  1.7201624
#> sigma:(Intercept)     sigma linpar (Intercept) -1.045460 0.07905694 -1.2004085
#>                        upper
#> mu:(Intercept)     1.2130200
#> mu:x               2.2955498
#> sigma:(Intercept) -0.8905109
confint(fit, "sigma")
#>                   parameter   term coefficient estimate         se     lower
#> sigma:(Intercept)     sigma linpar (Intercept) -1.04546 0.07905694 -1.200408
#>                        upper
#> sigma:(Intercept) -0.8905109
```
