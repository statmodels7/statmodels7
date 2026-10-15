# The Stacked Positions That a Specification Holds

Translates `spec@held_coef` into positions in the stacked coefficient
vector, from the design's own coefficient names, and adds the positions
that
[`term_held_stack()`](https://statmodels7.github.io/statmodels7/reference/term_held_stack.md)
returns for break-point terms held at fixed positions.

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
