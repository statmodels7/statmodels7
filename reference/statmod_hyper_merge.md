# Override the Starting Hyperparameters

Merges a caller's hyperparameters into the ones
[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
computed, by parameter and by term.

## Usage

``` r
statmod_hyper_merge(spec, start, user)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- start:

  The hyperparameters as
  [`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
  computed them.

- user:

  A named list of named lists, keyed by distribution parameter and then
  by term, or `NULL` for no overrides.

## Value

The same structure as `start`, with the caller's values merged in.
`start` unchanged when `user` is `NULL`.

## Details

Until an outer criterion estimates it, a hyperparameter sits at the
probe value, which is a placeholder: a lasso at \\\lambda = 1\\ against
an unaveraged log-likelihood over a few hundred observations selects
everything. This is the route by which a caller sets it instead.

A vector is matched by name against the penalty's own hyperparameters,
so `c(lambda = 5)` sets that one and leaves the rest. An unnamed vector
of the full length replaces them all.

A term may be named either by the key the specification holds it under,
which is its call deparsed, or by its `label`: `"lasso"` in place of
`"lasso(~noise1 + noise2)"`. Where two terms share a label the request
is ambiguous and the keys are asked for.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
for the structure,
[`hyper_key()`](https://statmodels7.github.io/statmodels7/reference/hyper_key.md)
for the name matching.
