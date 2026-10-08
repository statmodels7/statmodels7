# The Coordinates Left Out of the Marginal Determinant

The positions
[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md)
pins besides a boundary: the coefficients the criterion estimates, which
the specification holds, and the coordinates the model does not
identify, which
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
records on the design at its first fit.

## Usage

``` r
pinned_coords(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

An integer vector of stacked positions, possibly empty.
