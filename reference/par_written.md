# The Coordinates a Simulation's `par` Writes

Marks, for every distribution parameter and for the structural term's
own parameters, the positions that some key of `par` addresses.

## Usage

``` r
par_written(par, design, params, snm)
```

## Arguments

- par:

  A named list, or `NULL`.

- design:

  The design.

- params:

  The distribution parameter names.

- snm:

  The structural parameter names, or `NULL`.

## Value

A list with `beta`, one logical vector per distribution parameter, and
`zeta`, a logical vector over the structural parameters.

## Details

It reads the keys through
[`resolve_par_key()`](https://statmodels7.github.io/statmodels7/reference/resolve_par_key.md),
the function that
[`rstatmod_named()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_named.md)
uses to write them, so the two cannot disagree about what a key reaches.
A key that reaches nothing signals the same error here that it would
signal there.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md),
[`targets_written()`](https://statmodels7.github.io/statmodels7/reference/targets_written.md)
