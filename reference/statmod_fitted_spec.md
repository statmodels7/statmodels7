# The Terms as the Fit Left Them

Returns the specification with every refreshable term replaced by the
one the fit arrived at. A break-point, a nonlinear parameter and the
design block they imply are then read off the fitted object, and the
specification a caller passed to
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
no longer decides what
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) or
[`modelterms7::seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.html)
report.

## Usage

``` r
statmod_fitted_spec(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  the one the fit started from.

- coef:

  A named list of coefficient vectors, the ones the fit reached.

- design:

  The design at those coefficients, whose refresh state holds the terms
  to copy across.

## Value

A
[`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
identical to `spec` except that each refreshable term is the object the
fit left behind.

## Details

This is what a
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
stores in its `spec` property. The commit is `which = "jacobian"`: a
frozen block was committed by its own phase at exactly these
coefficients, and committing a
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
again at the same point would take a further step of its incremental
read-off, so the break-points reported would not be the fitted ones.
Both kinds of term are then copied across from the design's state,
committed or not.

A structural term's own parameters are copied too, from the design's
structural state onto `spec@structural`, so a filter's persistence and
loadings are read off the fit as well.

## See also

[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
for the commit this performs,
[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md)
for the terms replaced.
