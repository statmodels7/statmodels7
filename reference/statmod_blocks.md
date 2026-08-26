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
`sparse` (a list, one entry per non-smooth term, each with the
parameter, the term's name, its stacked columns and its penalty).

## Details

The property that decides is one each term already reports: a penalty
whose
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is non-empty is not twice differentiable in its coefficients, so its
block cannot enter a system solved by a curvature. Everything else goes
into one system and is estimated all together: an unpenalized block, a
ridge, a spline, a random effect, a structured or additive penalty.
Their joint curvature exists, and using it closes a fit in a handful of
iterations.

A term whose penalty gains a smooth approximation would answer
`penalty_kinks()` differently and move into the smooth block with no
change here.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
