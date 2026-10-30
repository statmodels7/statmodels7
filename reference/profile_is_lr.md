# Whether the Drop of a Profile Is a Likelihood Ratio

Twice the drop of the profile from its maximum is a likelihood ratio
statistic except in two cases: when a label is the hyperparameter of a
kinked penalty that a path chooses, whose profile includes the penalty
at the value held, and when the profile is the criterion
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and a label is a coefficient without a penalty, whose hold changes the
space of the fixed effects over which the criterion integrates.

## Usage

``` r
profile_is_lr(ctx, quantities, column)
```

## Arguments

- ctx:

  The context of the call.

- quantities:

  The
  [profiled_quantity](https://statmodels7.github.io/statmodels7/reference/profiled_quantity.md)
  objects of the labels.

- column:

  The column of the profile.

## Value

A single logical value.
