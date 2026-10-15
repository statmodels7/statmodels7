# Whether the Default Search Steers by the Exact Hessian

Whether a search this package chooses should take
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)'s
order-2 answer as its direction, which is a narrower question than
whether the Hessian exists.

## Usage

``` r
outer_newton_ok(spec, design, exact2, method = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, from
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md).

- exact2:

  Whether the criterion has an exact Hessian at all.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`. Where its `hessian` is `"expected"` the function returns
  `FALSE`.

## Value

`TRUE` or `FALSE`.

## Details

The two are separated because a criterion may have an exact second
derivative that is a poor direction for the search. Over a covariance
class the Hessian of the criterion is strongly indefinite wherever the
angle of the chart approaches its boundary, so that
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
repairs an indefinite Hessian and the repaired step can end elsewhere
than a quasi-Newton search. The default therefore uses
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
for a covariance class, whose Hessian is exact but not a reliable
direction.

The function returns `FALSE` for the expected information, for a model
with a structural term, and for a penalty with a covariance block wider
than one. In the first two cases the cost of the Hessian at each
evaluation outweighs the saving in evaluations.

The availability of the Hessian still governs the rest: it is supplied
to an optimizer that the caller names, it is what
[`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
and `vcov(type = "unconditional")` read, and it is recorded as
`exact_hessian`. Only the default direction is held back.

The block width is read from the `block` property of the penalty where
it has one; a penalty without that property is treated as univariate.

## See also

[`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md),
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
