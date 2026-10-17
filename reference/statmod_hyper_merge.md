# Override the Starting Hyperparameters

Merges a table of hyperparameter values, keyed by distribution parameter
and by term, into the structure that
[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
returns.

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
  returns them.

- user:

  A named list of named lists, keyed by distribution parameter and then
  by term, or `NULL` for no overrides.

## Value

The same structure as `start`, into which the caller's values are
merged. The result is `start` unchanged when `user` is `NULL`.

## Details

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does not call this function. It has no `hyper` argument, and a
hyperparameter is held in the term that carries the penalty, for example
`lasso(x, lambda = 3)`. The function is kept for code that builds a
hyperparameter structure by hand.

A vector is matched by name against the penalty's own hyperparameters,
so `c(lambda = 5)` sets that one and leaves the rest. An unnamed vector
of the full length replaces them all.

A term may be named either by the key under which the specification
holds it, which is its deparsed call, or by its `label`: `"lasso"` in
place of `"lasso(~noise1 + noise2)"`. Where two terms share a label the
name is ambiguous, and the function signals an error that lists the
keys.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
for the structure,
[`hyper_key()`](https://statmodels7.github.io/statmodels7/reference/hyper_key.md)
for the name matching.
