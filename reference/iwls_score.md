# The Dimensionless Reading of the Stopping Rule

The absolute \\\max\lvert g\rvert / n\\ read per equation against that
equation's own information scale, \\\max_j \lvert g_j\rvert / (n\\s_p)\\
with \\s_p = \sqrt{\mathrm{median}\_j H\_{jj} / n}\\ over the equation's
coordinates.

## Usage

``` r
iwls_score(g, hj, groups, n)
```

## Arguments

- g:

  The gradient at the point, a numeric vector over the stacked
  coefficients.

- hj:

  The information's diagonal at the same point, as
  [`iwls_info_diag()`](https://statmodels7.github.io/statmodels7/reference/iwls_info_diag.md)
  returns it, the same length as `g`.

- groups:

  A list of integer index vectors, one per equation, giving the
  coordinates that equation owns.

- n:

  The number of observations.

## Value

A single non-negative number: the largest per-equation reading. An
equation whose median curvature is zero or not finite contributes
nothing.

## Why the division is per equation

The score of a location equation carries the units \\1/y\\ and its
curvature \\1/y^2\\, so dividing by \\s_p\\ survives any rescaling of
the response. Taking \\s_p\\ over the whole equation changes nothing
within it: a stiff, heavily penalized coordinate is held to the same
rule as its neighbors, and that is what the envelope identities the
outer gradient rests on require of the mode.

This reading arbitrates the final verdict and never drives the loop.

## Two designs refused before this one

A per-coordinate normalization by \\(H+S)\_{jj}\\ is dimensionless too
and is wrong: at extreme shrinkage the penalized coordinates converge
loosely, the envelope identities lose their footing, and every outer
trajectory moved.

Driving the loop with the per-equation form is right at small \\y\\ and
unreachable at large. At \\y \times 10^4\\ the objective's magnitude
grows with \\\log y\\, so the stall guard on it fires before a rule
\\1/s_p\\ stricter can. The inner fit then reported failure across the
whole corridor of smoothing parameters between the criterion's plateau
and its optimum, and the outer search, which reads a non-converged point
as unavailable, never crossed it: 1482 evaluations and a fit at cor
0.82, where the absolute rule reaches 0.998.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
whose verdict this arbitrates,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the absolute rule that drives the loop.
