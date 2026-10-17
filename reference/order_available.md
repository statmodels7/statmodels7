# Whether a Family Carries Enough Derivatives

Returns `TRUE` when the log-density has at least `need` continuous
derivatives in every parameter for which the model carries an equation.
The consumers that can fall back on something else use the result as a
gate, and the consumers that cannot use it through
[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)
as an error message.

## Usage

``` r
order_available(distrib, need)
```

## Arguments

- distrib:

  A distributions7 family.

- need:

  The order required, an integer from 1 to 5: 1 for a score continuous
  in the coefficients, 2 for the curvature, 3 and 4 for the exact
  gradient and Hessian of a marginal criterion, and 5 for the fourth
  derivative of a filtered predictor.

## Value

A single logical. `NA` orders count as unavailable: a wrapper of a
kinked family records the kink without declaring its order, and an order
that has not been established is not evidence of smoothness.

## See also

[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)
for the message, and
[`distributions7::params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.html)
for the quantity.

## Examples

``` r
statmodels7:::order_available(distributions7::gaussian1_distrib(), 4L)
#> [1] TRUE
statmodels7:::order_available(distributions7::laplace_distrib(), 2L)
#> [1] FALSE
```
