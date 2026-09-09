# A Unit's Own Coordinates, in the Joint Vector's Order

The values whose positions
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
gives, returned in the order that function returns those positions in.

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
are a filter's own parameters and come from the design's structural
state, on the UNCONSTRAINED scale, which is the scale the recursion is
differentiated on and the one the penalty covers. A mixed class holds
some of each and interleaves them group by group, so its pieces are read
separately and put back in the class's own order, which is what
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md)
records.

Written beside
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
because a value and its address are one convention: a caller composing
them separately would scatter one unit's numbers at another's
coordinates, and the two would agree only by accident.

## See also

[`unit_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_beta.md),
the coefficient-only case this generalizes.
