# Refuse an Expected-Information Criterion That Would Read the Scores' Outer Product

Raises an error naming the family and the remedy when a criterion asked
for on the expected information would read, for this family, the outer
product of the scores at the data instead of an expectation.

## Usage

``` r
assert_criterion_information(distrib, method, approx)
```

## Arguments

- distrib:

  A distributions7 family.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`.

- approx:

  The approximation the fit reads the expected information with, from
  [`inner_settings()`](https://statmodels7.github.io/statmodels7/reference/inner_settings.md).

## Value

`NULL`, invisibly. Called for the error.

## Details

A criterion on the expected information replaces the observed curvature
of each observation,
\\-\partial^2\ell_i/\partial\eta\\\partial\eta^\top\\ at \\y_i\\, with
the Fisher information \\\mathcal{I}(\theta_i) =
-\mathbb{E}\[\partial^2\ell/\partial\eta\\ \partial\eta^\top\]\\, a
function of the parameters alone. Where a family writes that information
out,
[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
answers `TRUE` and nothing is asked here. Every shipped family does
since distributions7 0.65.0; where a family written outside the toolkit
does not, it falls back on the strategy `approx` names, and the default,
`"opg"`, returns \\s_i s_i^\top\\ with \\s_i\\ the score at \\y_i\\.
That is an unbiased estimate of \\\mathcal{I}(\theta_i)\\ from one
observation, not the information itself: it depends on the response, as
the observed curvature does, and it is neither of them. The criterion
built on it is neither the Laplace approximation nor its Fisher variant,
and it has no exact outer gradient, so its search is derivative-free and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
answers `unknown`.

Measured at 4000 observations with a smooth on the first parameter, on
the six families this applied to before distributions7 0.64.0 and 0.65.0
gave all of them an exact expected information, every one fit on the
observed information in 3 to 5 criterion evaluations with a Newton
search and an analytic certificate, where the outer-product route took
12 to 22 evaluations of a simplex.

The expectation is asked for by name with another `approx` in
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
and that is not refused: `"bartlett"` sums or integrates over the
support and is a genuine expected information, at a cost measured
between 66 and 1930 seconds per evaluation at 4000 observations.

## See also

[`expected_is_opg()`](https://statmodels7.github.io/statmodels7/reference/expected_is_opg.md),
the predicate, and
[`assert_criterion_order()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_order.md).

## Examples

``` r
statmodels7:::assert_criterion_information(
  distributions7::gaussian1_distrib(), reml("expected"), "opg")
# the Poisson-inverse Gaussian computes its expected information exactly
statmodels7:::assert_criterion_information(
  distributions7::pig1_distrib(), reml("expected"), "opg")
```
