# The Decrement Over the Directions Orthogonal to Some Others

The joint decrement over the coordinates `interior`, maximized over the
steps that have no component along the directions `dirs` in the
equilibrated coordinates.

## Usage

``` r
decrement_off(g, A, interior, dirs)
```

## Arguments

- g, A:

  The outer gradient and curvature.

- interior:

  The coordinates under test.

- dirs:

  A list of directions, vectors over every coordinate.

## Value

A single number, `NA` where the reduced curvature is not positive
definite.

## Details

A constraint on the step can only lower the maximum of \\2g'x - x'Ax\\,
so the result is at most
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
over the same coordinates, as removing a coordinate is.

## See also

[`certificate_ridges()`](https://statmodels7.github.io/statmodels7/reference/certificate_ridges.md),
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
