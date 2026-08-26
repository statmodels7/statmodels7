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

The fourth order is wanted where a filter's third derivative is: each
order of differentiating the predictor through the recursion pulls in
one more order of the family, the score the recursion is driven by being
read at the predictor it produces. Every family of distributions7
carries it in closed form.
