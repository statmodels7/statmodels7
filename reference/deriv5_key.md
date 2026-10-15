# The Name of a Fifth-Derivative Component

Locates the \\(a, b, c, d, e)\\ entry of a distribution's
fifth-derivative list, built the same way
[`deriv4_key()`](https://statmodels7.github.io/statmodels7/reference/deriv4_key.md)
builds its own.

## Usage

``` r
deriv5_key(params, a, b, c, d, e)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b, c, d, e:

  Indices into `params`.

## Value

A single string.

## Details

The fifth order is needed where a filter's fourth derivative is, by the
rule that the page of
[`deriv4_key()`](https://statmodels7.github.io/statmodels7/reference/deriv4_key.md)
states one order down: the criterion's second derivative in a pair of
hyperparameters reads
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
and the score that drives that recursion is read at the predictor that
the recursion produces.

[`distributions7::distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.html)
supplies the values in closed form for the families that register their
own fifth-order method, and as one central difference of the fourth
order for the others.
