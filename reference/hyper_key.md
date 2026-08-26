# Resolve a Term's Name Against a Specification

Turns the name a caller used into the key the specification holds the
term under, accepting either that key or the term's `label`.

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

  The distribution parameter whose equation to look in, a string.

- name:

  What the caller wrote: a key or a label.

## Value

A single string, the key `start[[p]]` holds the term under. Signals an
error when `name` matches nothing, and when it matches the labels of two
terms at once.

## Details

A specification keys its terms by the call as written, so a lasso is
`"lasso(~noise1 + noise2)"`. The call is what distinguishes two lassos
on different covariates, and it is also nothing anybody wants to type.
The `label` the term constructor carries is the short form, `"lasso"`,
and both are accepted here. Where two terms share a label the request is
ambiguous and the keys are asked for, guessing having a fair chance of
setting a hyperparameter on the wrong block.

## See also

[`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md),
its caller.
