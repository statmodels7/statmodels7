# Name Stacked Positions

The equation and the coefficient name of each position in the stacked
coefficient vector, in the shape
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
returns.

## Usage

``` r
stack_coords(spec, design, pos)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- pos:

  Integer positions in the stacked vector.

## Value

A list with `where`, `param` and `name`.
