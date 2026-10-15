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

A numeric vector of two, `lower` and `upper`. An end that the search
could not bracket within `maxit` widenings, or at which the statistic
could not be evaluated, is `NA`; both ends are `NA` when the standard
error of the coefficient is not finite and positive.

## Details

An interval by inversion is \\\\b : T(b) \le \chi^2\_{1,\alpha}\\\\, and
it is not symmetric about the estimate unless the statistic is a
parabola in \\b\\. The Wald statistic is a parabola by construction, the
curvature being read once at \\\hat\beta\\, and the other three
statistics are not. The asymmetry is the reason to invert one of the
other three: the interval follows the shape of the likelihood and not a
quadratic fitted at its maximum.

The search starts at the half-width of the Wald interval and widens the
distance from the estimate by a factor of 1.6, at most `maxit` times on
each side, until the statistic exceeds its critical value. It then
locates the crossing with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) (Brent's
method) at tolerance 1e-7. Each evaluation is one restricted refit,
warm-started at the unrestricted estimates and run at the fitted
hyperparameters, so the interval requires no outer search.

## See also

[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
the statistic being inverted.
