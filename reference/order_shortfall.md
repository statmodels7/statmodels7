# Why a Consumer Cannot Run on This Family

Builds the error message of a consumer that cannot run on the family.
The message names the family, the parameter, the order that the consumer
needs against the order that the family has, and what to write instead.
The result is `NULL` when the family has enough derivatives.

## Usage

``` r
order_shortfall(distrib, need, what)
```

## Arguments

- distrib:

  A distributions7 family.

- need:

  The order that the consumer needs.

- what:

  A noun phrase naming the consumer, used as the subject.

## Value

A single string, or `NULL`.

## Details

The remedy differs between two shapes of shortfall, and the message
distinguishes them. A family that DECLARES a kink has a known order and
the message names the composition that puts it there. A wrapper of such
a family declares none, so the order is not established rather than
known to be low, and the message reports that the order is not
established and gives no number.

## See also

[`order_available()`](https://statmodels7.github.io/statmodels7/reference/order_available.md),
[`distributions7::kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.html).

## Examples

``` r
statmodels7:::order_shortfall(distributions7::laplace_distrib(), 2L,
                              "the REML criterion")
#> [1] "the REML criterion needs 2 derivatives of the log-density in every parameter, and 'laplace' has 0 of them in 'mu' -- its log-density carries an absolute value of an argument that moves with 'mu', so the curvature that the REML criterion reads sits on a set of measure zero in the response. Fit with outer_criterion = NULL, which holds each smoothing hyperparameter at the value written in its term or at its starting value, or model 'mu' with a family that is smooth in it."
```
