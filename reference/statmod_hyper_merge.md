# Override the Starting Hyperparameters

Merges a caller's hyperparameters into the ones
[`statmod_hyper_start`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
computed, by parameter and by term.

## Usage

``` r
statmod_hyper_merge(spec, start, user)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- start:

  The hyperparameters as computed.

- user:

  A named list of named lists, or `NULL`.

## Value

The merged structure.

## Details

Until a hyperparameter is estimated by an outer criterion it is held at
the probe value, which is a placeholder and not a choice: a lasso at
\\\lambda = 1\\ against an unaveraged log-likelihood of a few hundred
observations selects nothing. This is what lets a caller set it.

A vector is matched by name against the penalty's own hyperparameters,
so `c(lambda = 5)` sets that one and leaves the rest where they were; an
unnamed vector of the full length replaces them all.

A term is named either by the key the specification holds it under,
which is its call deparsed, or by its `label` – `"lasso"` rather than
`"lasso(~noise1 + noise2)"`. Two terms sharing a label are ambiguous and
the keys are asked for instead.
