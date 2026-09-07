# A Confidence Interval by Inverting a Test

The set of values a test does not reject at `level`, found by locating
the two points where the statistic crosses its critical value.

## Usage

``` r
statmod_invert(
  fit,
  param,
  coefname,
  level = 0.95,
  test = "lr",
  type = "bayesian",
  maxit = 40L
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, coefname:

  Which coefficient, as in
  [`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md).

- level:

  The confidence level, a single number in (0, 1).

- test:

  One of `"lr"`, `"score"`, `"gradient"`, `"wald"`.

- type:

  Which variance sets the starting bracket, and which the Wald statistic
  reads.

- maxit:

  How many widenings to try on each side before giving up.

## Value

A numeric vector of two, `lower` and `upper`. An end the search could
not bracket is `NA`, which is what an unbounded interval looks like and
is reported rather than guessed at.

## Details

An interval by inversion is \\\\b : T(b) \le \chi^2\_{1,\alpha}\\\\, and
it is not symmetric about the estimate unless the statistic is a
parabola in \\b\\ – which the Wald statistic is BY CONSTRUCTION, the
curvature being read once at \\\hat\beta\\, and which the other three
are not. That asymmetry is the whole reason to invert one of the others:
it follows the likelihood's own shape rather than a quadratic fitted at
its top.

The search starts from the Wald interval, widens by a factor of 1.6
until the statistic exceeds its critical value, and then bisects. Each
evaluation is one restricted refit, warm-started at the unrestricted
estimates and run at the fitted hyperparameters, so an interval costs a
few dozen inner fits and no outer search at all.

## See also

[`statmod_stat()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat.md),
the statistic being inverted.
