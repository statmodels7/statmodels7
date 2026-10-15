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

  Which coefficients: a distribution parameter's name, a term's name, a
  vector of `parameter:coefficient` labels, or `NULL` for all of them.
  With `readable = TRUE`, the name of a quantity is accepted as well.

- level:

  The confidence level.

- type:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  `"bayesian"`, `"frequentist"` or `"unconditional"`.

- readable:

  Whether to report the quantities in which the model is written (a
  term's own parameters under the names of the literature) and not the
  raw coordinates. Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
  and necessarily `FALSE` for the three restricted tests.

- test:

  Which test the interval is the acceptance region of: `"wald"`, the
  default, or `"lr"`, `"score"` or `"gradient"`, each inverted by
  [`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md).

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A data frame with the parameter, the term, the coefficient, the
estimate, its standard error and the two limits.

## Details

With `readable = FALSE`, the interval for a coefficient of a linear
predictor is symmetric about the estimate and needs no mapping back.
Such a coefficient is unbounded whatever the distribution parameter it
belongs to, the link having already carried that parameter onto the
whole line, so the scale on which the interval is built is the scale on
which the quantity lives. The interval does not respect a bound on the
distribution parameter itself; for that, map an interval for the
predictor through the inverse link at the covariate values of interest.
With the default `readable = TRUE`, each quantity that a structural term
declares (a loading or a persistence, for example) has an interval built
on the link scale of that quantity and mapped back, so that interval is
not symmetric.

The variance comes from
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
so the same three conventions apply, and a coefficient that a kinked
penalty has set to zero has `NA` in place of an interval. For an
interval around a penalized term, `type = "unconditional"` does not
treat the smoothing parameter as known.

## Inverting a test instead

`test` chooses which statistic the interval is the acceptance region of.
The Wald interval is the set of \\b\\ that the Wald statistic does not
reject, and it is symmetric by construction: that statistic is a
parabola in \\b\\, the curvature having been read once at \\\hat\beta\\.
The other three read the likelihood at each \\b\\ in turn, so their
intervals follow its shape and are not symmetric where it is not. Where
the likelihood is nearly quadratic the four intervals agree; with a
small sample and a skewed likelihood the likelihood-ratio and gradient
intervals shift away from the Wald interval.

Each of the three is found by
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which brackets the crossing outward from the Wald interval and locates
it with a root finder, so one interval requires a few dozen restricted
refits. The `parm` argument therefore selects the coefficients before
the search and not after it, so that only the requested intervals are
computed.

A row for which the restricted fit is not defined keeps `NA` limits.
[`testable_coords()`](https://statmodels7.github.io/statmodels7/reference/testable_coords.md)
identifies those rows: a coefficient under a kinked penalty, an aliased
one, and every coefficient of a model carrying a structural term. The
estimate and the standard error in such a row are those of the fit and
do not depend on `test`.

A coefficient under a penalty that is twice differentiable (a smooth's
own coordinates, a ridge, a random effect) is inverted like any other,
and the interval is the credible one that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
gives for it under `type = "bayesian"`, conditional on the
hyperparameters that the fit reached. See
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).

The three restricted tests (`"lr"`, `"score"` and `"gradient"`) need
`readable = FALSE`, since a test is about one coefficient and a readable
quantity is a function of several at once.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
which tests one coefficient against a value of its own, and
[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md)
and
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
the two internal functions from which an inverted interval is built

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
confint(fit)
#>                   parameter   term coefficient  estimate         se      lower
#> mu:(Intercept)           mu linpar (Intercept)  1.043291 0.08770104  0.8714004
#> mu:x                     mu linpar           x  2.007856 0.14865515  1.7164973
#> sigma:(Intercept)     sigma linpar (Intercept) -1.032801 0.08006408 -1.1897235
#>                        upper
#> mu:(Intercept)     1.2151822
#> mu:x               2.2992148
#> sigma:(Intercept) -0.8758781
confint(fit, "sigma")
#>                   parameter   term coefficient  estimate         se     lower
#> sigma:(Intercept)     sigma linpar (Intercept) -1.032801 0.08006408 -1.189723
#>                        upper
#> sigma:(Intercept) -0.8758781
confint(fit, "mu:x", test = "lr", readable = FALSE)
#>      parameter   term coefficient estimate        se    lower    upper
#> mu:x        mu linpar           x 2.007856 0.1486552 1.716674 2.299038
```
