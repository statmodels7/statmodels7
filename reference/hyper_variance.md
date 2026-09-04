# The Variance of the Estimated Hyperparameters

Inverts the negative of the outer criterion's Hessian, holding any
coordinate whose own curvature cannot produce a variance.

## Usage

``` r
hyper_variance(A, schur = 1e-04)
```

## Arguments

- A:

  The negative of the outer Hessian, with dimnames.

- schur:

  The largest relative Schur correction a held coordinate may contribute
  to a kept one's curvature.

## Value

A matrix of the same shape as `A` with the variance in the kept rows and
columns and `NA` elsewhere, or `NULL` when no coordinate is usable or
the coupling is too large to ignore.

## Details

At a maximum the criterion's Hessian is negative definite and its
negative inverts to a variance. A hyperparameter driven to the edge of
its range, or one the search left before reaching a maximum, has a
curvature there that is zero or of the wrong sign, and no variance
follows from it.

Such a coordinate used to cost every other one its standard error, the
whole matrix being refused. It is held instead and the rest is inverted,
which is the variance CONDITIONAL on it — the same reading
[`vcov()`](https://rdrr.io/r/stats/vcov.html) gives a coefficient the
information carries nothing about. That is the marginal variance only
where the coupling contributes nothing to the kept curvature, so the
Schur correction \\A\_{kb}A\_{bb}^{-1}A\_{bk}\\ is computed and compared
against the kept diagonal rather than assumed negligible; above `schur`
the whole matrix is refused as before. The default is the size at which
the correction cannot move the four significant digits the summary
prints.

## See also

[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md),
its only caller.
