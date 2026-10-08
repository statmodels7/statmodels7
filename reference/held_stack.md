# The Stacked Positions a Specification Holds

Translates `spec@held_coef` into positions in the stacked coefficient
vector, from the design's own coefficient names.

## Usage

``` r
held_stack(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

An integer vector of positions, empty where nothing is held.
