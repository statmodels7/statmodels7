# The Distribution Parameters That a Marginal Criterion Estimates

Resolves the `marginal` argument of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
against the family: the parameters whose unpenalized coefficients are
estimated by maximizing the marginal criterion.

## Usage

``` r
marginal_params(method, distrib)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`.

- distrib:

  The family.

## Value

A list with `named`, a character vector of parameter names in the
family's order, and `explicit`, a single logical that is `TRUE` when the
caller named them and `FALSE` when the default did. When the caller
named parameters the list also has `all`, `TRUE` for `marginal = "all"`.

## Details

`NULL` names every parameter except the position, which is the family's
first parameter. In most shipped univariate families the first parameter
is the location (the `mu` of
[`distributions7::weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.html)).
A few families start with another kind of parameter: the first parameter
of
[`distributions7::beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.html)
is the shape `alpha` and that of
[`distributions7::gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.html)
is the scale `sigma`.

A multivariate family is not covered. An explicit request signals an
error, and the default names nothing, which is the convention of the
criterion without the argument.

## See also

[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
which turns the names into coefficients.
