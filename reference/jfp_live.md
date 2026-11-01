# Small Helpers of the Third and Fourth Derivatives

`jfp_live()` returns the equations in which a block has rows;
`jfp_add_weights()` sums weight lists with coefficients, and
`jfp_add_pair()` adds a trace term of
[`jfp_mv_pair()`](https://statmodels7.github.io/statmodels7/reference/jfp_mv_ctx.md)
with a coefficient; `jfp_darray()` returns the family's derivative of
the expected information of order `k`, kept in the state;
`jfp_check_outer()` signals an error where a block's information is not
positive definite.

## Usage

``` r
jfp_live(bs)

jfp_add_weights(acc, lst, coef)

jfp_add_pair(acc, pair, coef)

jfp_darray(st, k)

jfp_check_outer(st)
```

## Arguments

- bs:

  A block's entry of the state.

- acc, lst, coef:

  An accumulated weight list, the lists to add and their coefficients.

- pair:

  What
  [`jfp_mv_pair()`](https://statmodels7.github.io/statmodels7/reference/jfp_mv_ctx.md)
  returns.

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- k:

  The order, 1 to 4.

## Value

As described.
