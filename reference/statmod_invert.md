# A Confidence Interval by Inverting a Test

Returns the set of values that a test does not reject at `level`, found
by locating the two points where the statistic crosses its critical
value.

## Usage

``` r
statmod_invert(
  fit,
  param,
  coefname,
  level = 0.95,
  test = lr_test(),
  vcov_type = bayesian_vcov(),
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

  The
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
  to invert:
  [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md)
  (the default),
  [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md),
  [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)
  or
  [`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md).

- vcov_type:

  The
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  that sets the starting bracket, and that the Wald statistic reads.

- maxit:

  How many widenings to try on each side before the search stops.

## Value

A numeric vector of two, `lower` and `upper`. A side of the position of
a break-point that the data do not bound is `-Inf` or `Inf`. An end that
the search could not bracket within `maxit` widenings, or at which the
statistic could not be evaluated, is `NA`; both ends are `NA` when the
standard error of the coefficient is not finite and positive.

## Details

An interval by inversion is \\\\b : T(b) \le \chi^2\_{1,\alpha}\\\\, and
it is not symmetric about the estimate unless the statistic is a
parabola in \\b\\. The Wald statistic is a parabola by construction,
because the curvature is read once at \\\hat\beta\\, and the other three
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

The position of a break-point has a range: a value at or beyond the
smallest or the largest observation of the variable leaves the term with
no effect, so that every such value gives the same refit.
[`break_edges()`](https://statmodels7.github.io/statmodels7/reference/break_edges.md)
returns that range, and the search never leaves it. On a side where the
statistic is still under its critical value when the search reaches the
end of the range, the data do not bound the position, and the limit is
`-Inf` or `Inf`.

## See also

[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md),
the statistic that is inverted.
