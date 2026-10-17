# A Unit's Own Coordinates, in the Joint Vector's Order

Returns the values whose positions
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
gives, in the order in which that function returns those positions.

## Usage

``` r
unit_joint_beta(u, spec, design, coef)
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

- coef:

  The coefficients, one entry per distribution parameter.

## Value

A numeric vector as long as
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)'s
answer.

## Details

The three kinds are read from two different places. An ordinary unit's
coordinates are coefficients and come from `coef`. A structural unit's
coordinates are a filter's own parameters and come from the design's
structural state, on the unconstrained scale, which is the scale on
which the recursion is differentiated and which the penalty covers. A
mixed class holds some of each and interleaves them group by group, so
its pieces are read separately and put back in the class's own order,
which
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md)
records.

This function is written beside
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
because a value and its position must follow the same convention. A
caller that composed them separately could place one unit's values at
another unit's positions.

## See also

[`unit_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_beta.md),
the coefficient-only case that this function generalizes.
