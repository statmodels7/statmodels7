# The Name of a Third-Derivative Component

Locates the \\(a, b, c)\\ entry of a distribution's third-derivative
list, which is keyed by
[`distributions7::deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.html).

## Usage

``` r
deriv3_key(params, a, b, c)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b, c:

  Indices into `params`.

## Value

A single string.

## Details

The name is built from the parameter names in the family's own order and
is not parsed out of an existing name, so that a parameter whose own
name contains the separator is handled as distributions7 handles it.
