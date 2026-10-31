# The Derivatives of the Expected Information Kept in a Jeffreys State

`jfp_dexpected()` returns
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and `jfp_d2expected()` returns
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html),
both on the link scale at the state's predictors, computed at the first
request and kept in the state.

## Usage

``` r
jfp_dexpected(st)

jfp_d2expected(st)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

## Value

A named list of per-observation components.
