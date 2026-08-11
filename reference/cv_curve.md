# The Held-Out Deviance of Every Point of a Path

Refits the model on each training fold along the whole path and scores
it on the fold left out, returning the mean deviance per observation and
its standard error across folds.

## Usage

``` r
cv_curve(spec, data, weights, offsets, inner_method, hypers, folds)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- data:

  The data the fit was called on.

- weights, offsets:

  As
  [`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  received them.

- inner_method:

  The inner method.

- hypers:

  A list of hyperparameter settings, one per path point.

- folds:

  A fold number per observation.

## Value

A list with `cvm`, `cvse` and `n_fail`.

## Details

The path is run fold by fold rather than point by point so that each fit
starts from the previous point's coefficients, which is what makes a
path cheaper than its length suggests. Each training fit rebuilds the
design on its own rows: a term is re-evaluated in the data it is fitted
to, so a basis or a set of contrasts is not carried over from rows the
fit did not see.

## See also

[`cv`](https://statmodels7.github.io/statmodels7/reference/cv.md)
