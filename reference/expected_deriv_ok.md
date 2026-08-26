# Does the Family Supply the Expected Information's Derivative?

Whether
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
answers for this family, asked at a probe and never inferred from its
class.

## Usage

``` r
expected_deriv_ok(distrib)
```

## Arguments

- distrib:

  A distributions7 distribution.

## Value

A single logical.

## Details

The default route in distributions7 is one central difference of the
family's own expected information, which is a single stencil on an
analytic quantity wherever that information is a written-out formula and
refuses where it is itself an integral. Six of forty univariate families
refuse, and the reason is cost, never accuracy: measured at 100
observations they cost 1880 to 147300 ms against a median of 0.183 ms
for the others, so 2p of those calls per criterion evaluation is not a
slower route but an unusable one.

## See also

[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
