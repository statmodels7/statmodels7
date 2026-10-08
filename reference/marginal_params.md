# The Distribution Parameters a Marginal Criterion Estimates

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
family's order, and `explicit`, a single logical saying whether the
caller named them (`TRUE`) or the default did (`FALSE`).

## Details

`NULL` names every parameter except the position, which is the family's
first parameter. In every shipped univariate family the first parameter
is the location or, where the family has none, the scale (the `mu` of
[`distributions7::weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.html),
the `sigma` of
[`distributions7::gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.html)).

A multivariate family is not covered. An explicit request is refused and
the default names nothing, which is the convention of the criterion
without the argument.

## See also

[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
which turns the names into coefficients.
