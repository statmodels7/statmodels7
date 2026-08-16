# Evaluate the Offsets a Formula Names

The offset per distribution parameter, evaluated in the data.

## Usage

``` r
eval_offsets(formula, params, data, env, n)
```

## Arguments

- formula:

  The model formula, before the offsets are stripped.

- params:

  The distribution's parameter names.

- data:

  A data frame.

- env:

  The environment the formula carried.

- n:

  The number of observations.

## Value

A named list, one entry per parameter, `NULL` where the equation names
no offset.

## Details

The expressions are re-evaluated rather than carried as numbers, which
is what lets an offset survive prediction:
[`statmod_respec`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
calls this against the new data, where a vector supplied through the
`offsets` argument at fitting time has the wrong length and cannot be
reused.

## See also

[`split_offsets`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
