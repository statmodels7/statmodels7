# Does the Family's Expected Information Depend on the Response?

The result is `TRUE` where
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

The answer comes from a probe and not from a list of families, so that a
family written later is covered. The probe draws two distinct responses
at one parameter value and compares the two answers. A family that
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
reports as writing its information out is not probed, and the probe is
not run for any `approx` other than `"opg"`, the only strategy that
reads the response. For a truncated family, the predicate returns
`FALSE` although the expected information is a quadrature that does not
depend on the response, so this function also returns `FALSE`, and the
answer is determined by the probe and not by the predicate. The caller's
random stream is restored.

## See also

[`assert_criterion_information()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_information.md)
