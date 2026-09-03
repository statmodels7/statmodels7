# The Joint Prior of a Covariance Class

The penalty over a class's stacked coefficients: the distribution one of
its members named, or the centered multivariate Gaussian that is the
default.

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

The prior describes the effects of **one group** over every column the
label collects, so its dimension is the class's total and the penalty
covers \\m\\ such blocks. Which chart the covariance rides is a modeling
choice and is
[`parameters7::dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.html)
by default, where a coordinate is the logarithm of a standard deviation
exactly; the alternative, log-Cholesky, is the same family of matrices
written so that only the first coordinate reads as one.

Naming a `distrib` on more than one member is an error rather than a
precedence rule: the prior is one object and there is nothing to say
which of two should win.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md),
its only caller.
