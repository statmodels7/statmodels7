# One Penalized Unit, by Parameter and Key

Returns the entry of
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
that a hyperparameter row names, or `NULL` where there is none.

## Usage

``` r
statmod_unit(spec, design, param, key)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- param:

  The distribution parameter.

- key:

  The key, as
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  composes it.

## Value

One entry, or `NULL`.

## Details

A reader that starts from a `(parameter, term)` pair looks the penalty
up here and not from the term's whole block, which assumes one penalty
per term. The lookup returns that penalty where the term has one, and
the penalty with the given key where it has several.
