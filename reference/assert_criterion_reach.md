# Prediction-Error Criteria and a Structural Term's Own Parameters

Raises an error naming the penalty and the remedy when
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
would have to select a hyperparameter of a penalty over the own
parameters of a structural term: the group departures of a score-driven
filter over a panel, or a covariance class that reaches into one.

## Usage

``` r
assert_criterion_reach(
  spec,
  design,
  method,
  reach = c("all", "smooth", "kinked")
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`.

- reach:

  Which penalties the criterion selects: `"all"`, `"smooth"` or
  `"kinked"`.

## Value

`NULL`, invisibly. Called for the error.

## Details

A prediction-error criterion reads the effective degrees of freedom, and
their derivative, over the coefficients alone, where such a penalty
covers nothing. The exact gradient of
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
there is zero, so a search would stop at its first evaluation and report
convergence.
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) and
the path along which a kinked penalty is swept are built on the same
blocks.
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
span the term's parameters in their determinant and estimate a twice
differentiable penalty there. A penalty with a kink there has no
criterion that can select its hyperparameter. The hyperparameter can be
held by writing its value in the term, for example
`lasso(~id, lambda = 2)`, and the penalty is then fitted over the term's
own parameters by a proximal route, without a criterion.

Which penalties a criterion reaches depends on its role, and only those
are checked: an `outer_criterion` reaches the smooth penalties, and the
kinked ones as well when no `sparse_criterion` is given, while a
`sparse_criterion` reaches the kinked ones alone.

## See also

[`assert_criterion_order()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_order.md),
[`structural_penalized()`](https://statmodels7.github.io/statmodels7/reference/structural_penalized.md)

## Examples

``` r
dd <- data.frame(x = runif(60))
dd$y <- rnorm(60, dd$x)
spec <- statmodels7:::statmod_spec(y ~ x, distributions7::gaussian1_distrib(),
                                   dd)
statmodels7:::assert_criterion_reach(spec, statmodels7:::statmod_design(spec),
                                     aic())
```
