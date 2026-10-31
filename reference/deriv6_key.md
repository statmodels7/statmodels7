# The Name of a Sixth-Derivative Component

Locates the \\(a, b, c, d, e, f)\\ entry of a distribution's
sixth-derivative list, built the same way
[`deriv5_key()`](https://statmodels7.github.io/statmodels7/reference/deriv5_key.md)
builds its own.

## Usage

``` r
deriv6_key(params, a, b, c, d, e, f)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b, c, d, e, f:

  Indices into `params`.

## Value

A single string.

## Details

The sixth order is read where a filter's fifth derivative is,
[`modelterms7::term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.html),
for the reason that
[`deriv5_key()`](https://statmodels7.github.io/statmodels7/reference/deriv5_key.md)
gives one order down.
[`distributions7::distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.html)
supplies the values.
