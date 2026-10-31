# The Jeffreys State at Given Coefficients, Kept Between Calls

Returns
[`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md)
for the Jeffreys blocks of the model at `coef`, or `NULL` where the
model carries no Jeffreys term. The objective, its gradient and its
curvature read the same point in turn, so the last state is kept on the
design's memo.

## Usage

``` r
jfp_state_at(spec, design, coef, keys = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- keys:

  As in
  [`jfp_blocks()`](https://statmodels7.github.io/statmodels7/reference/jfp_blocks.md);
  a state over some of the units is not kept.

## Value

An environment, or `NULL`.
