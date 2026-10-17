# Where a Penalized Unit's Coordinates Are in the Joint Vector

Returns the same positions as
[`unit_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_positions.md),
carried onto the one vector \\\[\beta; \zeta\_{\mathrm{free}}\]\\ that
the inner step, the marginal criterion and the variance all assemble, so
that a unit of any of the three kinds can be read off one diagonal.

## Usage

``` r
unit_joint_positions(u, spec, design)
```

## Arguments

- u:

  One unit, from
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

An integer vector of positions, empty where the unit names a coordinate
that the joint vector does not carry.

## Details

The three kinds are addressed in three different vectors, and a caller
that forms a trace needs one. A mixed class already records the answer,
because
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md)
builds its `joint` field. A structural unit records positions among all
of a term's parameters, where the joint vector holds only the free ones,
so a held coordinate has to be mapped out rather than offset past. An
ordinary unit's stacked index is already a position in the joint vector,
because the coefficients come first and every one of them is a free
coordinate.

This function is written once so that the convention is composed in one
place and not in three.

## See also

[`unit_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_positions.md)
for the same coordinates in the unit's own vector;
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md),
which builds the mixed case.
