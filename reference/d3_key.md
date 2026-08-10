# The Name of a Third-Derivative Component

Locates the \\(a, b, k)\\ entry of a distribution's third derivative,
which is keyed by name and not by position.

## Usage

``` r
d3_key(params, a, b, k, keys)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b, k:

  Indices into `params`.

- keys:

  The names the derivative actually returned.

## Value

A single string.

## Details

The name is BUILT by putting the three parameter names in the family's
own order and joining them, the direction distributions7 sanctions, and
then checked against the enumeration rather than trusted.
