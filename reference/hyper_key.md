# Resolve a Term's Name Against a Specification

Returns the key under which the specification holds a term, given either
that key or the `label` of the term.
[`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md)
is the only caller, and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does not use it.

## Usage

``` r
hyper_key(spec, start, p, name)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  read for the terms' labels.

- start:

  The hyperparameter structure, whose names are the keys.

- p:

  The distribution parameter in whose equation to look, a string.

- name:

  What the caller wrote: a key or a label.

## Value

A single string, the key under which `start[[p]]` holds the term. The
function signals an error when `name` matches nothing, and when it
matches the labels of two terms at once.

## Details

A specification keys its terms by the call as written, so a lasso is
`"lasso(~noise1 + noise2)"`. The call distinguishes two lassos on
different covariates. The `label` of a term is a shorter name,
`"lasso"`, and both are accepted. Where two terms share a label the name
is ambiguous, and the function signals an error that lists the keys.

## See also

[`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md),
its caller.
