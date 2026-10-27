# The Family's Second-Derivative Components at a Point

Computes the Hessian of the distribution, or its expected information,
on the link scale at one value of the parameters.

## Usage

``` r
statmod_family_hessian(spec, theta, expected = TRUE, emethod = NULL)
```

## Arguments

- spec:

  The specification, read for the distribution, the response and the
  thread count.

- theta:

  The parameters, one value per observation.

- expected:

  `TRUE` for the expected information, `FALSE` for the observed Hessian.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

The named list of second-derivative components on the link scale, as
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html)
returns it when `expected` is `TRUE` and as
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
otherwise.

## Details

The components are computed in one place because several routes read
them at the same point.
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
builds a square-root decomposition from the per-observation blocks and
uses the assembled information when that decomposition cannot be formed.
Both routes receive the components computed here, so the family is
evaluated once per point.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md),
[`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md),
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
and
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md),
which read it.
