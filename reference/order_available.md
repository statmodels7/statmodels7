# Whether a Family Carries Enough Derivatives

`TRUE` when the log-density has at least `need` continuous derivatives
in every parameter the model carries an equation for. Used as a gate by
the consumers that can fall back on something else, and through
[`order_shortfall()`](https://statmodels7.github.io/statmodels7/reference/order_shortfall.md)
as a refusal by the ones that cannot.

## Usage

``` r
order_available(distrib, need)
```

## Arguments

- distrib:

  A distributions7 family.

- need:

  The order asked for; see the table in `R/order.R`.

## Value

A single logical. `NA` orders count as unavailable: a wrapper of a
kinked family records the kink without declaring its order, and an
unestablished order is not evidence of smoothness.

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
