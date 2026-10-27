# Settle the Curvature of an `iwls()` Against a Family

Turns the information of an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
into the matrix that a fit runs. With `information` given, that object
replaces the one that the method holds; without it, a method that holds
an unsettled
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
is settled and a method whose curvature was named is returned as it is.

## Usage

``` r
iwls_resolve(method, distrib, information = NULL)
```

## Arguments

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- distrib:

  The distribution to which the model is fitted.

- information:

  `NULL`, or an object that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html)
  that replaces the one of `method`.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object whose `hessian` is an
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
or an
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
object, unless a named choice is returned as it is.

## Details

A named
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
or
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
is the matrix of the step and nothing else: the method carries no
fallback and no switch.

An
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
is settled by two questions.
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is the first. Where it returns `FALSE` the expected information is an
approximation, the outer product of the scores by default, and the
method is settled on the observed information, with the expected
information as the `backup` that takes the step where the observed step
finds no acceptable point. An observed penalized information that is not
positive definite is first floored by
[`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
and its step tried.

[`distributions7::expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.html)
is the second question. A family whose expected information is exact but
costly is settled on the observed information as well, because the
expected information of such a family requires a quadrature or a sum
over the support at each observation. The expected information is the
`backup` there too.

Only a family settled on the expected information has `switching` set.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which reads `backup`.
