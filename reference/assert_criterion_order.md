# Refuse Where a Criterion Cannot Be Computed

Raises
[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)'s
message when the family does not carry the two derivatives every
criterion reads. Called from
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
once the criterion is known to have something to do, so a model whose
criterion is inert is not refused for a reason that would never have
arisen.

## Usage

``` r
assert_criterion_order(distrib, method)
```

## Arguments

- distrib:

  A distributions7 family.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`.

## Value

`NULL`, invisibly. Called for the error.

## See also

[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)

## Examples

``` r
statmodels7:::assert_criterion_order(distributions7::gaussian1_distrib(), reml())
```
