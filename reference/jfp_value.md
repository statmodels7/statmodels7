# The Value of the Jeffreys Penalty

The Value of the Jeffreys Penalty

## Usage

``` r
jfp_value(st)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

## Value

\\\sum_b a_b\log\|\tilde I_b\|\\, a single number, `-Inf` where a
block's information is not positive definite.

## See also

[`jfp_gradient()`](https://statmodels7.github.io/statmodels7/reference/jfp_gradient.md),
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
