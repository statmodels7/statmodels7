# Does the Family's Expected Information Depend on the Response?

`TRUE` where
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html),
called with `approx`, returns a matrix that changes with the response at
fixed parameters, which an expectation cannot do.

## Usage

``` r
expected_is_opg(distrib, approx)
```

## Arguments

- distrib:

  A distributions7 family.

- approx:

  The approximation, a string.

## Value

A single logical.

## Details

Asked at a probe rather than read from a list of families, so a family
written later is covered: two distinct responses drawn at one parameter
value, and the two answers compared. A family that
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
reports as writing its information out is never probed, and neither is
any `approx` other than `"opg"`, the only strategy that reads the
response. A truncated family reports `FALSE` to that predicate although
its expected information is a quadrature that does not depend on the
response, so it answers `FALSE` here: that is why the probe, and not the
predicate, decides. The caller's random stream is restored.

## See also

[`assert_criterion_information()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_information.md)
