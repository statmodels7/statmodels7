# Does the Family Supply the Expected Information's Second Derivative?

Whether
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html)
answers for this family, asked at a probe.

## Usage

``` r
expected_deriv2_ok(distrib)
```

## Arguments

- distrib:

  A distributions7 distribution.

## Value

A single logical.

## Details

There is no numerical default in distributions7 – a difference of the
first derivative, itself a difference for most families, would be the
nested differencing that package forbids – so the families that answer
are the ones that wrote the derivative out: `gaussian1_distrib()`,
`poisson_distrib()`, `gamma1_distrib()`, `negbin2_distrib()` and
`beta1_distrib()`. Every other family keeps the expected route's order 2
on
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md).

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md),
[`expected_deriv_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv_ok.md)
