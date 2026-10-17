# The Held-Out Criterion of Every Point of a Path

Refits the model on each training fold along the whole path and scores
it on the fold left out. Returns, at each point, the mean over folds of
minus twice the weighted log-likelihood of the held-out rows divided by
the sum of their weights, and its standard error across folds.

## Usage

``` r
cv_curve(spec, data, inner_optimizer, hypers, folds, run = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  holding the weights and offsets of the fit.

- data:

  The data with which the fit was called.

- inner_optimizer:

  The inner method.

- hypers:

  A list of hyperparameter settings, one per path point.

- folds:

  A fold number per observation.

- run:

  The combination of the outer axes to which each point belongs. The
  warm start begins again at the head of each combination, where the
  kink jumps back up. `NULL` treats the whole list as one run.

## Value

A list with `cvm`, the mean over folds at each point; `cvse`, its
standard error across the folds with a finite value; `n_fail`, the
number of folds without a value at each point; and `error`, the message
of the first fold that could not rebuild the model, or `NULL`. `cvm` is
missing at a point where every fold lacks a value.

## Details

The path is run fold by fold and not point by point, so that each fit
starts from the previous point's coefficients. Each fold reapplies to
its rows the terms built on all the rows, as
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
reapplies them to new data: the knots of a basis, the levels of a factor
and the scales of a standardized block come from the covariates of every
row and never from the response. A basis rebuilt on the training rows
alone would not cover a held-out row beyond their range, and the fold
would fail.

The training rows and the held-out rows carry the weights and offsets of
`spec` at those rows, including the offsets given through the `offsets`
argument of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).
Where the fit of a fold at a point does not converge, or its criterion
is not finite, the fold's value at that point is missing and is counted
in `n_fail`.

## See also

[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
