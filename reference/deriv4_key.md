# The Name of a Fourth-Derivative Component

Locates the \\(a, b, c, d)\\ entry of a distribution's fourth-derivative
list, built the same way
[`deriv3_key()`](https://statmodels7.github.io/statmodels7/reference/deriv3_key.md)
builds its own.

## Usage

``` r
deriv4_key(params, a, b, c, d)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b, c, d:

  Indices into `params`.

## Value

A single string.

## Details

The fourth order is needed where a filter's third derivative is: each
order of differentiating the predictor through the recursion requires
one more order of the family, because the score that drives the
recursion is read at the predictor that the recursion produces. The
values come from
[`distributions7::distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.html).
