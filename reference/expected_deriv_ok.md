# Does the Family Supply the Expected Information's Derivative?

Determines whether
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
returns a value for this family. The generic is called once, at a probe
value, and the result is not inferred from the family's class.

## Usage

``` r
expected_deriv_ok(distrib, emethod = NULL)
```

## Arguments

- distrib:

  A distributions7 distribution.

## Value

A single logical.

## Details

The default method in distributions7 is one central difference of the
family's expected information. It applies where that information is a
written-out formula and signals an error where the information is itself
approximated, because differencing an approximation would be a
difference of a difference. The function returns `FALSE` for a family on
which the generic signals an error, and also where a random parameter
vector or a response cannot be drawn for the probe.

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
