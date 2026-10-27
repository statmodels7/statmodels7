# The Information of a Variance Convention for a Family

Returns the `information` property of a
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
object that has one and that is not `NULL`, and
[`distributions7::default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.html)
of `distrib` otherwise (a `NULL` property, or a convention that has no
property).

## Usage

``` r
vcov_information(x, distrib)
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

## Value

An object that inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## See also

[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
[`distributions7::default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.html)
