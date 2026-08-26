# Evaluate the Offsets a Formula Names

The offset per distribution parameter, evaluated in the data.

## Usage

``` r
eval_offsets(formula, params, data, env, n)
```

## Arguments

- formula:

  The model formula, before
  [`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
  has stripped anything.

- params:

  The distribution's parameter names, in the family's order.

- data:

  A data frame to evaluate in.

- env:

  The environment the formula carried, the enclosure of the evaluation.

- n:

  The number of observations, the length to recycle to.

## Value

A named list with one entry per element of `params`, each a numeric
vector of length `n` or `NULL` where that equation names no offset.

## Details

The expressions are kept and re-evaluated, never carried as numbers, and
that is how an offset survives prediction.
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
calls this against the new data; a vector supplied through
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`offsets` argument at fitting time has the wrong length for other rows
and cannot be reused. Before this, `predict(fit, newdata =)` returned
the predictor of a model with no offset at all.

Each expression is evaluated in `data` with `env` behind it, so a symbol
resolves as a column first and as a variable of the caller's environment
second. A result shorter than `n` is recycled, so a single number is a
constant offset.

## See also

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
