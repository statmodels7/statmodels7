# The Variance of the Hyperparameters a Marginal Criterion Estimated

The asymptotic variance matrix of the estimated hyperparameters, on the
free scale their links carry them onto.

## Usage

``` r
statmod_hyper_vcov(
  spec,
  design,
  coef,
  hyper,
  method,
  inner = NULL,
  pinned = NULL,
  at_edge = character(0)
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  The outer method that estimated them, or `NULL`.

- inner:

  The inner optimizer that the fit used, with which the stencil route
  refits its probes;
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  where none is given.

- pinned:

  The coordinates that the fit left out of the criterion's determinant,
  as it records them in `methods$pinned`.

- at_edge:

  The `boundary_key` of
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md):
  the hyperparameters at the edge of their range, held by
  [`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)
  so that the others carry the variance conditional on them.

## Value

A square matrix, one row per estimated hyperparameter, whose dimnames
join the distribution parameter, the term and the hyperparameter's own
name with a carriage return, a character that no name can contain, so
that the three are recoverable from the key. The index is stored in the
attribute `"idx"`. The rows and columns of held coordinates are `NA`
(see
[`hyper_variance()`](https://statmodels7.github.io/statmodels7/reference/hyper_variance.md)).
The value is `NULL` where there is nothing to report.

## Details

A hyperparameter estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
is the maximizer of a criterion that is twice differentiable in it, so
its variance is the inverse of the negative Hessian of that criterion at
the point reached, which
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
returns. If the criterion also estimates coefficients, the variance is
the block of the hyperparameters in the inverse of the negative Hessian
over the hyperparameters and those coefficients.
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
computes the Hessian analytically where a route applies and by a
difference of the exact gradient otherwise. The matrix is read on the
free scale because the criterion was maximized on that scale, where the
quadratic approximation is reasonable. A variance on the scale of the
hyperparameter itself would describe a shape that the criterion does not
have where the estimate is far from zero and the criterion is
asymmetric.

**Where it does not apply.** The function returns `NULL` for every outer
criterion other than
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
including
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
used as the outer criterion for a smooth penalty. A hyperparameter that
a path selects (the penalty of a kinked term, chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) over
a grid) is the argument of a minimum, not a root of a derivative, so no
Hessian exists to invert. Its uncertainty can be assessed by resampling.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
