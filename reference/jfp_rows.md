# The Rows of a Jeffreys Block

Returns \\\tilde V_q\\, one matrix per distribution parameter, `NULL`
where the block does not enter that parameter's predictor.

## Usage

``` r
jfp_rows(st, b)
```

## Arguments

- st:

  The state under construction, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- b:

  A
  [`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md).

## Value

A list as long as the distribution's parameters.

## Details

For a block of kind `"static"` the matrix of its equation is the block's
columns of the design, preceded by a column of ones when the constant
direction is included. A sparse design is densified over the block's
columns only.
