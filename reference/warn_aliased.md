# Report a Model Whose Columns Are Not All Identified

Raises a warning naming the coordinates the fit could not estimate, so
that a reader who does not print the summary is told the model as
written does not identify every column of its design.

## Usage

``` r
warn_aliased(alias)
```

## Arguments

- alias:

  The aliased coordinates by name, as
  [`aliased_labels()`](https://statmodels7.github.io/statmodels7/reference/aliased_labels.md)
  gives them.

## Value

`NULL`, invisibly. Called for the warning.

## Details

A coefficient the model does not identify is one point of a flat ridge:
the fit converges and reports numbers, and nothing in what it prints by
default says which of them mean anything.
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) put `NA` in
[`coef()`](https://rdrr.io/r/stats/coef.html) and say no more, and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
has followed that convention since 0.98.0, but an `NA` among a hundred
rows of a [`summary()`](https://rdrr.io/r/base/summary.html) nobody
printed is not a report.

Which coordinates those are is
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)'s
answer, read at the fitted coefficients on \\K = H + S\\ and gated by
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
and that is what makes the warning safe to raise. A rank test on the RAW
design answers a different question and has false positives that would
make this noise: measured over fourteen models it reports a deficiency
on `random(~ 1 | g)`, on `random(~ x | g)` and on
`s(x) + random(~ 1 | g)`, whose columns are identified by their own
penalty and whose fits are perfectly ordinary. Read on \\K\\, the
penalty's curvature is in the matrix and none of the three is named.
Over the same fourteen the warning fires on exactly two, both of them
`nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))`, where a parameter's own
intercept sits beside a complete set of indicators that sum to it.

Reading it at the FITTED coefficients rather than at the start is part
of the rule: a Jacobian block is degenerate at a starting value of zero
for a reason that goes away as soon as the parameter moves, and
`nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))`, which is identified, is
deficient there and is not named here.

THE MESSAGE NAMES NO CAUSE, and a measurement is why. Two shapes reach
it and they are not the same defect. In the one above the DESIGN is
deficient, a column being the sum of others. In the other the design has
full rank and \\K\\ does not: measured on
`jump(x, smoothed = smooth_quintic())` and on `jseg()` under the same
smoother, the block is of rank 3 of 3 while `mu:jump.psi1` is named, the
quintic being exact outside its own width so that the break-point's
derivative column is non-zero at a handful of observations and carries
almost no curvature. Both are genuine – neither coefficient has a
standard error – and a message asserting the first would be wrong about
the second, so it reports the fact and leaves the diagnosis to the
reader. The two smoothers whose derivative has unbounded support,
`smooth_probit()` and `smooth_hyperbolic()`, are not named on the same
data.

## See also

[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md),
which answers which coordinates those are, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which reports them as missing.
