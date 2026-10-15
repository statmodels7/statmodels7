# Does the Family Supply the Expected Information's Second Derivative?

Whether
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html)
returns a value for this family. The generic is called once, at a probe
value.

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

The default method in distributions7 signals an error, because
differencing the first derivative of the expected information would be a
difference of a difference. The function therefore returns `TRUE` only
for a family that provides the derivative itself. For any other family
the order-2 Hessian of the expected route is computed by
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md).

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md),
[`expected_deriv_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv_ok.md)
