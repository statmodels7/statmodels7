# The Variance of the Hyperparameters a Marginal Criterion Estimated

The asymptotic variance matrix of the estimated hyperparameters, on the
free scale their links carry them onto.

## Usage

``` r
statmod_hyper_vcov(spec, design, coef, hyper, method)
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

## Value

A square matrix, one row per estimated hyperparameter, whose dimnames
join the distribution parameter, the term and the hyperparameter's own
name with a carriage return, a character no name can contain, so the
three stay recoverable from the key. The index rides on the attribute
`"idx"`. `NULL` where there is nothing to report.

## Details

A hyperparameter estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
is the maximizer of a criterion that is twice differentiable in it, so
it has a variance like any other maximum-likelihood estimate: the
inverse of the negative Hessian of that criterion at the point reached,
which
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
already computes exactly. It is read on the free scale because that is
where the criterion was maximized and where the quadratic approximation
behind it is reasonable; a variance for a smoothing parameter on its own
scale, where the estimate is often several orders of magnitude from zero
and the criterion far from symmetric, would describe a shape the
criterion does not have.

**Where it does not apply.** A hyperparameter chosen along a path, which
is what
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
do to a kinked penalty's, is the argument of a minimum over a grid. It
is no root of a derivative, so there is no Hessian to invert and no
standard error follows. Its uncertainty is a resampling question, and
`NULL` is returned in place of a number of another kind.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
