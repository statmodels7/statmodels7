# The Constant Direction and the Power of Each Jeffreys Unit

Completes the units of
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
that carry a
[`modelterms7::JeffreysPenalty()`](https://statmodels7.github.io/modelterms7/reference/JeffreysPenalty.html):
each receives `power`, the power of its determinant, and `const`,
whether the constant direction of its equation enters the determinant.

## Usage

``` r
jfp_resolve_units(units, spec, design)
```

## Arguments

- units:

  The units, as
  [`statmod_penalized_units()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  builds them.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

`units`, with `power` and `const` set on each Jeffreys unit.

## Details

The constant direction follows the rule of
[`jfp_constant()`](https://statmodels7.github.io/statmodels7/reference/jfp_constant.md)
with the setting `intercept` of the term. The level's columns are the
columns of the equation that no penalized unit covers and that do not
move with the coefficients: the columns of a term that recomputes its
block (the Jacobian of
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
or of a smoothed
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html))
are left out, so that the decision is the same at every point at which
the specification is built. With only such a term in the level the
constant direction is never included.

## See also

[`jfp_blocks()`](https://statmodels7.github.io/statmodels7/reference/jfp_blocks.md)
