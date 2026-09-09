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

The fifth order is wanted where a filter's FOURTH derivative is, by the
rule the fourth's own page states one order down: the criterion's second
derivative in a pair of hyperparameters reads
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
and the score that recursion is driven by is read at the predictor it
produces.

It is the one place the fifth order enters, and
[`distributions7::distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.html)
supplies it as ONE central difference of the analytic fourth rather than
in closed form. What that costs is measured on the page of
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md).
