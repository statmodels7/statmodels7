# Split a Specification's Terms Into the Smooth Block and the Rest

Returns which coefficients are fitted jointly and which carry a penalty
with a kink and are therefore fitted by a method of their own.

## Usage

``` r
statmod_blocks(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A list with `smooth` (the stacked column indices fitted jointly) and
`sparse` (a list with one entry per penalized unit that has a kink).
Each entry holds `param`, `term` (the penalty's key), `cols` (positions
within the parameter's coefficients), `index` (positions in the stacked
vector), `penalty`, `structural`, `zterm` (the term's own name) and
`ids` (the labels under which hyperparameters are shared).

## Details

The split depends on the penalty of each penalized unit (see
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)).
A penalty for which
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
returns a finite point is not twice differentiable in its coefficients,
so its block cannot enter a system solved by a curvature. Everything
else goes into one system and is estimated all together: an unpenalized
block, a ridge, a spline, a random effect, a structured or additive
penalty. Their joint curvature exists, so a scoring or Newton step
applies to them.

The classification is made by
[`penalty_has_kink()`](https://statmodels7.github.io/statmodels7/reference/penalty_has_kink.md),
so a penalty that becomes smooth moves into the smooth block with no
change here.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
