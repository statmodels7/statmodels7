# The Jeffreys Blocks of a Model

Returns one
[`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
per Jeffreys unit of the model, in the order of
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

## Usage

``` r
jfp_blocks(spec, design, keys = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- keys:

  `NULL` for every Jeffreys unit, or the keys of the units to keep.

## Value

A list of
[`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
objects, empty where the model carries no Jeffreys term.
