# The Coordinates That a Kinked Penalty Holds at Zero

The positions among those of
[`kinked_coords()`](https://statmodels7.github.io/statmodels7/reference/kinked_coords.md)
whose coefficient equals zero.

## Usage

``` r
zero_kinked(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- coef:

  The coefficients, a named list.

## Value

An integer vector of stacked positions, possibly empty.
