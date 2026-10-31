# Small Helpers of the Third and Fourth Derivatives

`jfp_live()` returns the equations in which a block has rows;
`jfp_add_weights()` sums weight lists with coefficients; `jfp_darray()`
returns the family's derivative of the expected information of order
`k`, kept in the state; `jfp_check_outer()` signals an error where the
third and fourth derivatives are not written.

## Usage

``` r
jfp_live(bs)

jfp_add_weights(acc, lst, coef)

jfp_darray(st, k)

jfp_check_outer(st)
```

## Arguments

- bs:

  A block's entry of the state.

- acc, lst, coef:

  An accumulated weight list, the lists to add and their coefficients.

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- k:

  The order, 1 to 4.

## Value

As described.
