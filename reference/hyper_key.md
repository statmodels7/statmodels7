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
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- start:

  The hyperparameter structure, whose names are the keys.

- p:

  The distribution parameter.

- name:

  What the caller wrote.

## Value

A single key.

## Details

A specification keys its terms by the call as written, so a lasso is
`"lasso(~noise1 + noise2)"`. That is what makes two lassos on different
covariates distinct, and it is not what anybody wants to type; the label
the term constructor carries is. Where two terms share a label the
request is ambiguous and the keys are asked for, since guessing would
set a hyperparameter on the wrong block.
