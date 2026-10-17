# The Joint Prior of a Covariance Class

Returns the penalty over a class's stacked coefficients: the
distribution that one of its members named, or the centered multivariate
Gaussian that is the default.

## Usage

``` r
class_penalty(cl)
```

## Arguments

- cl:

  One class, as
  [`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md)
  assembles it before the penalty.

## Value

A penalties7 penalty over `m * dim` coefficients.

## Details

The prior describes the effects of one group over every column that the
label collects, so its dimension is the class's total and the penalty
covers \\m\\ such blocks. The default prior is a centered multivariate
Gaussian whose covariance is written by
[`parameters7::dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.html),
where a coordinate is the logarithm of a standard deviation exactly.

Naming a `distrib` on more than one member is an error and not a
precedence rule: the prior is one object, so two of them cannot both be
used.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md),
its only caller.
