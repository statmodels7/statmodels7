# Settle the Curvature of an `iwls()` Against a Family

Turns an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
left at `hessian = "auto"` into the method a fit runs: the expected
information where the family writes it exactly and at an acceptable
cost, and otherwise the observed information, with the expected one
taking the step where the observed step finds no acceptable point. A
method whose curvature was named is returned as it is.

## Usage

``` r
iwls_resolve(method, distrib)
```

## Arguments

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- distrib:

  The distribution the fit is of.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object whose `hessian` is `"expected"` or `"observed"`, with `fallback`
`TRUE` exactly where `"auto"` was settled on the observed information.

## Details

[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is the first question. Where it answers `FALSE` the expected information
is an approximation, the outer product of the scores by default, and the
method is settled on the observed information. An observed penalized
information that is not positive definite is first repaired by
[`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
and its step tried; the expected information takes the step only where
that step finds no acceptable point.

[`distributions7::expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.html)
is the second question. A family whose expected information is exact but
costly is settled on the observed information as well, because the
expected information of such a family requires a quadrature or a sum
over the support at each observation. The expected information is still
used where the observed one cannot take the step, and in full where a
method or a criterion names it.

Only a family settled on the expected information receives a
`switch_after` other than `Inf`; see
[`iwls_switch_after()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_after.md).

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which reads `fallback`.
