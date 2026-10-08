# The Held-Out Deviance of Every Point of a Path

Refits the model on each training fold along the whole path and scores
it on the fold left out, returning the mean deviance per observation and
its standard error across folds.

## Usage

``` r
cv_curve(
  spec,
  data,
  weights,
  offsets,
  inner_optimizer,
  hypers,
  folds,
  run = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- data:

  The data the fit was called on.

- weights, offsets:

  As
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  received them.

- inner_optimizer:

  The inner method.

- hypers:

  A list of hyperparameter settings, one per path point.

- folds:

  A fold number per observation.

- run:

  Which combination of the outer axes each point belongs to. The warm
  start begins again at the head of each, the kink jumping back up
  there. `NULL` treats the whole list as one run.

## Value

A list with `cvm`, `cvse` and `n_fail`.

## Details

The path is run fold by fold, not point by point, so that each fit
starts from the previous point's coefficients. That warm chain is the
whole economy of a path cheaper than its length suggests. Each fold
reapplies to its rows the terms built on all the rows, as
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
reapplies them to new data: the knots of a basis, the levels of a factor
and the scales of a standardized block come from the covariates of every
row, never from the response. Rebuilt on the training rows alone, a
basis did not cover a held out row past their range, and the fold
stopped the whole fit (a smooth of the year of construction beside a
lasso on `gamlss.data::rent`).

## See also

[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
