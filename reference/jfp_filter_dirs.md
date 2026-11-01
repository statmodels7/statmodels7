# The Directions of a Jeffreys Block Inside a Filter

Returns the block's directions over every coordinate of the joint vector
(the coefficients followed by all the filter's parameters, held levels
included): the columns of `C` placed at the block's positions among the
filter's parameters, so that the block's rows are the derivatives of the
predictors along these directions.

## Usage

``` r
jfp_filter_dirs(st, b)
```

## Arguments

- st:

  The state under construction, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- b:

  A
  [`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
  of kind `"filter"`.

## Value

A numeric matrix, the coordinates by the block's columns.
