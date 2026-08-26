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

The name is built from the parameter names in the family's own order,
not parsed out of one, which is the discipline distributions7 records
for a parameter whose own name contains the separator.
