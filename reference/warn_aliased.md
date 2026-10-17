# Report a Model Whose Columns Are Not All Identified

Raises a warning naming the coordinates that the fit could not estimate,
so that a reader who does not print the summary learns that the model as
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

`NULL`, invisibly. The function is called for the warning.

## Details

A coefficient that the model does not identify is one point of a flat
ridge: the fit converges and reports numbers, and the default print does
not indicate which of them are meaningful.
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) put `NA` in
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
follows that convention. The warning also reports the coordinates when
[`summary()`](https://rdrr.io/r/base/summary.html) is not printed, where
an `NA` among the rows of a table is easy to overlook.

The coordinates are those of
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md),
read at the fitted coefficients on \\K = H + S\\ and gated by
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
which makes the warning safe to raise. A rank test on the raw design
addresses a different question and has false positives that would make
the warning noisy: it reports a deficiency on `random(~ 1 | g)`, on
`random(~ x | g)` and on `s(x) + random(~ 1 | g)`, whose columns are
identified by their own penalty and whose fits are ordinary. Read on
\\K\\, the curvature of the penalty is in the matrix and the function
names none of the three. The warning does fire on
`nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))`, where the intercept of a
parameter sits beside a complete set of indicators that sum to it.

Reading the matrix at the fitted coefficients and not at the start is
part of the rule: a Jacobian block is degenerate at a starting value of
zero for a reason that goes away as soon as the parameter moves, and
`nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))`, which is identified, is
deficient there and is not named here.

The message names no cause, because two different shapes reach it. In
the first the design is deficient because a column is the sum of others.
In the second the design has full rank and \\K\\ does not, as for
`jump(x, smoothed = smooth_quintic(h = g))` at one median spacing \\g\\,
and for `jseg()` under the same smoother: the block has rank 3 of 3
while `mu:jump.psi1` is named, because the quintic is exact outside its
own width, so that where every observation lies farther than \\h\\ from
the break-point its derivative column carries no curvature. The quintic
resolves its width at \\5/(2\log 2)\\ median spacings, which avoids this
on a uniform covariate; a width supplied by hand can still produce it.
Both cases are genuine (in each of them the coefficient has no standard
error), and a message asserting the first would be wrong about the
second, so the message reports the fact and leaves the diagnosis to the
reader. The two smoothers whose derivative has unbounded support,
`smooth_probit()` and `smooth_hyperbolic()`, are not named on the same
data.

## See also

[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md),
which determines which coordinates those are, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which reports them as missing.
