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
equation whose median curvature is zero or not finite is scaled by one,
so that its score per observation enters the reading as it is.

## Why the division is per equation

The score of a location equation carries the units \\1/y\\ and its
curvature \\1/y^2\\, so dividing by \\s_p\\ survives any rescaling of
the response. Taking \\s_p\\ over the whole equation changes nothing
within it: a stiff, heavily penalized coordinate is held to the same
rule as its neighbors, and the envelope identities on which the outer
gradient rests require this of the mode.

This reading arbitrates the final verdict and never drives the loop.

## Alternative readings that were rejected

A per-coordinate normalization by \\(H+S)\_{jj}\\ is dimensionless too,
but at extreme shrinkage the penalized coordinates converge loosely and
the envelope identities lose their footing.

Using the per-equation form to drive the loop works at small \\y\\ and
cannot be met at large \\y\\. At \\y \times 10^4\\ the magnitude of the
objective grows with \\\log y\\, so the stall guard on the objective
stops the run before a rule \\1/s_p\\ stricter can be met.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
whose verdict this arbitrates,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the absolute rule that drives the loop.
