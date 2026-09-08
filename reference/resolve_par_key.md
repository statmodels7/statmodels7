# Resolve One Key of a Simulation's `par`

Turns a name into the vector and the positions it addresses.

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

The order is exact matches first and groups after, so a name that is
both a coefficient and the stem of others addresses the coefficient. A
group is recognized at a dot rather than by any prefix, since a
coefficient name is free to begin with the letters of another.

A key that reaches nothing is reported with what the model does carry.
Guessing would be worse than refusing: a misspelled name would leave the
quantity drawn and the caller would read a simulation of another model.

## See also

[`rstatmod_named()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_named.md)
