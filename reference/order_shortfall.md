# Why a Consumer Cannot Run on This Family

The message a consumer refuses with: the family, the parameter, the
order it needs against the order there is, and what to write instead.
`NULL` when there is nothing to refuse.

## Usage

``` r
order_shortfall(distrib, need, what)
```

## Arguments

- distrib:

  A distributions7 family.

- need:

  The order the consumer needs.

- what:

  A noun phrase naming the consumer, used as the subject.

## Value

A single string, or `NULL`.

## Details

Two shapes of shortfall are told apart, because the remedies differ. A
family that DECLARES a kink has a known order and the message says which
composition puts it there. A wrapper of such a family declares none, so
the order is unestablished rather than known to be low, and the message
says that instead of inventing a number.

## See also

[`order_available()`](https://statmodels7.github.io/statmodels7/reference/order_available.md),
[`distributions7::kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.html).

## Examples

``` r
statmodels7:::order_shortfall(distributions7::laplace_distrib(), 2L,
                              "the REML criterion")
#> [1] "the REML criterion needs 2 derivatives of the log-density in every parameter, and 'laplace' has 0 of them in 'mu' -- its log-density carries an absolute value of an argument that moves with 'mu', so the curvature the REML criterion reads sits on a set of measure zero in the response. Fit with outer_criterion = NULL, setting the hyperparameters yourself through hyper, or model 'mu' with a family that is smooth in it."
```
