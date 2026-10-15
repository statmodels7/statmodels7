# Resolve One Key of a Simulation's `par`

Turns a name into the vector and the positions that it addresses.

## Usage

``` r
resolve_par_key(key, params, cn, snm)
```

## Arguments

- key:

  The name.

- params:

  The distribution parameter names.

- cn:

  A named list of coefficient names, one per parameter.

- snm:

  The structural parameter names, or `NULL`.

## Value

A list with `space`, `param` and `pos`.

## Details

A key is tried in this order: a distribution parameter; a structural
parameter, then a group of structural parameters; a coefficient of an
equation, written `parameter.coefficient`, then a group of the
coefficients of that equation. A name that is both a coefficient and the
stem of others therefore addresses the coefficient. A group is
recognized at a dot and not by any prefix, since a coefficient name may
begin with the letters of another.

A key that addresses nothing is reported with what the model does carry.
A misspelled name that was ignored would leave the quantity drawn at
random, and the simulation would be of another model than the one
intended.

## See also

[`rstatmod_named()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_named.md)
