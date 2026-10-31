# The Method of the Expected Information That a Jeffreys Penalty Reads

Returns `NULL` where the family has the closed form of its expected
information and of its derivatives, and
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
where it has only its density.

## Usage

``` r
jfp_emethod(distrib)
```

## Arguments

- distrib:

  A
  [distributions7::distrib](https://statmodels7.github.io/distributions7/reference/distrib.html)
  object.

## Value

`NULL` or an
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
object.

## Details

The penalty is the Jeffreys prior only on the expected information, so
the method of the fit is not read. The outer product of the scores at
the data, which
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
computes, is not the expected information, and a family that reaches it
is rejected; so is a family on which
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
gives no first and second derivative of the matrix.
