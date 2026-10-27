# The Information of a Variance Convention for a Family

Returns the `information` property of a
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
object that has one and that is not `NULL`, and
[`distributions7::default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.html)
of `distrib` otherwise (a `NULL` property, or a convention that has no
property).

## Usage

``` r
vcov_information(x, distrib, design = NULL)
```

## Arguments

- x:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- distrib:

  The
  [distributions7::distrib](https://statmodels7.github.io/distributions7/reference/distrib.html)
  of the fit.

- design:

  The design of the fit, or `NULL` where it is not at hand.

## Value

An object that inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## Details

Where `design` carries a structural term of the likelihood shape
([`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html),
a marginal break-point term), the default is
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
whatever the family, and a property that names
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
signals an error (see
[`structural_expected_check()`](https://statmodels7.github.io/statmodels7/reference/structural_expected_check.md)):
the likelihood of such a model is a mixture over latent components and
has no expected information.

## See also

[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
[`distributions7::default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.html)
