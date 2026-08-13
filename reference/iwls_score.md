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

  The gradient at the point.

- hj:

  The information's diagonal, from
  [`iwls_info_diag`](https://statmodels7.github.io/statmodels7/reference/iwls_info_diag.md).

- groups:

  The equations' coordinate index sets.

- n:

  The number of observations.

## Value

A single number.

## Details

The score of a location equation carries the units \\1/y\\ and its
curvature \\1/y^2\\, so the division survives any rescaling of the
response, while changing nothing WITHIN an equation: a stiff, heavily
penalized coordinate is held to the same rule as its neighbours, which
is what the envelope identities the outer gradient rests on ask of the
mode. It arbitrates the final verdict only. Two designs were tried and
refused before this one: a per-coordinate normalization by
\\(H+S)\_{jj}\\ let the penalized coordinates converge loosely at
extreme shrinkage and moved every outer trajectory, and driving the LOOP
with the per-equation form made the tolerance unreachable at \\y \cdot
10^4\\, where the stall guard on the objective, whose magnitude grows
with \\\log y\\, fires before a rule \\1/s_p\\ stricter can – the inner
then reported failure across the whole corridor of smoothing parameters
between the plateau and the optimum, and the outer search, reading those
points as unavailable, never crossed it (1482 evaluations, a fit at cor
0.82 where the absolute rule's run reaches 0.998).
