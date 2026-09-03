# The Coefficients a Penalized Unit Covers

The values the unit's penalty is read at, in the unit's own order.

## Usage

``` r
unit_beta(u, coef, params)
```

## Arguments

- u:

  One unit, from
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- coef:

  The coefficients, a named list by distribution parameter.

- params:

  The distribution's parameters, in order.

## Value

A numeric vector as long as the unit's index.

## Details

Written once because a unit's coefficients are addressed differently
depending on what it is, and every caller wants the same vector. An
ordinary unit sits in one equation and its columns are positions in that
parameter's coefficients; a covariance class spans several and is
addressed only in the stacked vector. Stacking and indexing answers
both, and for an ordinary unit it returns exactly what
`coef[[u$param]][u$cols]` returned, the stacked index being that
parameter's offset plus those columns.

A structural unit has no position in the stacked vector at all: its
penalty covers the term's own parameters, which contribute no design
column. It is read from the design's structural state instead and never
reaches here.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md).
