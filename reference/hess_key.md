# The Name of a Second-Derivative Component

Locates the \\(a, b)\\ entry of a distribution's Hessian list, which is
keyed by
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.html)
and not by position.

## Usage

``` r
hess_key(params, a, b)
```

## Arguments

- params:

  The parameter names, in the family's order.

- a, b:

  Indices into `params`.

## Value

A single string.

## Details

The name is BUILT from the parameter names rather than parsed out of
one, the discipline distributions7 records for a parameter whose own
name contains the separator.
