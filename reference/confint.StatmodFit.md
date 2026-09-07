# Confidence Intervals for a Fit

Confidence intervals for the coefficients of every distribution
parameter, by default the Wald ones and otherwise by inverting a
likelihood test.

## Usage

``` r
# S3 method for class 'StatmodFit'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("bayesian", "frequentist", "unconditional"),
  readable = TRUE,
  test = c("wald", "lr", "score", "gradient"),
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
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  `"bayesian"`, `"frequentist"` or `"unconditional"`.

- readable:

  Whether to report the quantities the model is written in – a term's
  own parameters under the names its literature gives them – rather than
  the raw coordinates. Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
  and necessarily `FALSE` for the three restricted `method`s.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

- method:

  Which test the interval is the acceptance region of: `"wald"`, the
  default, or `"lr"`, `"score"` or `"gradient"`, each inverted by
  [`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md).

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
so the same three conventions apply, and a coefficient a kinked penalty
set to zero has `NA` in place of an interval. An interval around a
penalized term is worth asking for as `type = "unconditional"`, which
does not read the smoothing parameter as though it had been known.

## Inverting a test instead

`method` chooses which statistic the interval is the acceptance region
of. The Wald interval is the set of \\b\\ the Wald statistic does not
reject, and it is symmetric BY CONSTRUCTION: that statistic is a
parabola in \\b\\, the curvature having been read once at \\\hat\beta\\.
The other three read the likelihood at each \\b\\ in turn, so their
intervals follow its own shape and are not symmetric where it is not.
Measured on a Poisson regression at \\n = 500\\, where the likelihood is
nearly quadratic, the four agree to three decimals; at \\n = 25\\ with a
skewed fit the likelihood-ratio interval sits \\0.046\\ further right
than the Wald one and the gradient interval \\0.071\\.

Each of the three is found by
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which brackets from the Wald interval and bisects, so one interval is a
few dozen restricted refits – measured, about six times one fit of the
model. It is therefore asked for one coefficient at a time in practice,
through `parm`, and the subsetting happens BEFORE the search rather than
after so that nothing is paid for that is not reported.

A row the restricted fit is not defined for keeps `NA` limits:
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
says which those are – a coefficient under a KINKED penalty, an aliased
one, and every coefficient of a model carrying a structural term. The
estimate and the standard error in such a row are the fit's own and
stand whatever `method` is.

A coefficient under a penalty that is twice differentiable – a smooth's
own coordinates, a ridge, a random effect – is inverted like any other,
and the interval is the credible one
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
already gives it under `type = "bayesian"`, conditional on the
hyperparameters the fit reached. See
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).

The three restricted methods need `readable = FALSE`, since a test is
about one coefficient and a readable quantity is a function of several
at once.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
which tests one coefficient against a value of its own, and
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md)
and
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
the two internals an inverted interval is built from

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
confint(fit, "mu:x", test = "lr", readable = FALSE)
#>      parameter   term coefficient estimate        se    lower    upper
#> mu:x        mu linpar           x 2.007856 0.1467852 1.716674 2.299038
```
