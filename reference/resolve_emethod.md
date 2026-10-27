# The Method That Obtains the Derivatives of the Expected Information

Returns the method object with which the derivatives of the expected
information (the functions
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html))
are called, so that the derivative is the derivative of the matrix that
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html)
returns for the same `emethod`.

## Usage

``` r
resolve_emethod(distrib, emethod)
```

## Arguments

- distrib:

  A
  [distributions7::distrib](https://statmodels7.github.io/distributions7/reference/distrib.html)
  object.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

`NULL` or an object that inherits from
[distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).

## Details

An object is returned unchanged. For `NULL`, the matrix is the closed
form of the family where the family has one, and the derivatives of the
family are those of that closed form, so `NULL` is returned; where the
family has no closed form the matrix is
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
and the derivatives signal an error with `NULL`, so
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
is returned and used for both. A family is read as having the closed
form of its derivative when
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is `TRUE`, and otherwise (a truncated family, whose expected information
is a quadrature that the family registers with its derivative) when
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
returns a value at a random point; the probe leaves the random stream of
the caller as it was and the answer for the last family probed is kept.

## See also

[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md),
[`expected_deriv_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv_ok.md)
