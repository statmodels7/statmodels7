# One Penalized Unit, by Parameter and Key

The entry of
[`statmod_penalized`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
a hyperparameter row names, or `NULL` where there is none.

## Usage

``` r
statmod_unit(spec, design, param, key)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

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

The places that read a penalty from a `(parameter, term)` pair used to
fetch it with `term_penalty()` and take the term's whole block, which
assumes one penalty per term. Looking it up here answers the same
question where that holds and the right question where it does not.
