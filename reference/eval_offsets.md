# Evaluate the Offsets That a Formula Names

Returns the offset of each distribution parameter, evaluated in the
data.

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

  The environment that the formula carried, the enclosure of the
  evaluation.

- n:

  The number of observations.

## Value

A named list with one entry per element of `params`, each a numeric
vector of length `n` or `NULL` where that equation names no offset.

## Details

The expressions are kept and re-evaluated, not carried as numbers, so
that an offset written in the formula applies in prediction.
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
calls this against the new data; a vector supplied through
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`offsets` argument at fitting time has the length of the fitting data
and cannot be reused on other rows.

Each expression is evaluated in `data` with `env` behind it, so a symbol
resolves as a column first and as a variable of the caller's environment
second. A result of length one is recycled to `n`, so a single number is
a constant offset, and a result of any other length except `n` is an
error. A value that is missing or not finite is an error. Several
offsets in one equation are summed.

## See also

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
