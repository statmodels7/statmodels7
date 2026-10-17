# The Coefficients That a Penalized Unit Covers

Returns the values at which the unit's penalty is read, in the unit's
own order.

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

This function is written once because a unit's coefficients are
addressed differently depending on what the unit is, and every caller
needs the same vector. An ordinary unit sits in one equation and its
columns are positions in that parameter's coefficients; a covariance
class may span several equations and is addressed only in the stacked
vector. Stacking and indexing serves both, and for an ordinary unit it
returns exactly what `coef[[u$param]][u$cols]` returned, because the
stacked index is that parameter's offset plus those columns.

A structural unit has no position in the stacked vector at all: its
penalty covers the term's own parameters, which contribute no design
column. It is read from the design's structural state instead. A mixed
covariance class is addressed in the joint vector.
[`unit_joint_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_beta.md)
reads both of them.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md).
