# Derivative Orders Read by a Criterion

Raises
[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)'s
message when the family does not carry the two derivatives that every
criterion reads.
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
calls the function once the criterion is known to have something to do,
so that a model whose criterion is inert does not signal an error that
cannot arise.

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

`NULL`, invisibly. The function is called for the error.

## See also

[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)

## Examples

``` r
statmodels7:::assert_criterion_order(distributions7::gaussian1_distrib(), reml())
```
