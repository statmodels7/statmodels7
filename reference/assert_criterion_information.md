# Expected-Information Criteria and the Scores' Outer Product

Raises an error naming the family and the remedy when a criterion
defined on the expected information would read, for this family, the
outer product of the scores at the data instead of an expectation.

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
answers `TRUE` and no check is made. Where a family does not write it
out, for example a family written outside the toolkit, the information
falls back on the strategy `approx` names, and the default, `"opg"`,
returns \\s_i s_i^\top\\ with \\s_i\\ the score at \\y_i\\. That is an
unbiased estimate of \\\mathcal{I}(\theta_i)\\ from one observation, not
the information itself: it depends on the response, as the observed
curvature does, and it is neither of them. The criterion built on it is
neither the Laplace approximation nor its Fisher variant, and it has no
exact outer gradient, so its search is derivative-free and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
answers `unknown`.

A criterion on the observed information has an exact outer gradient and
Hessian and is searched with Newton's method; the outer-product route
has neither and is searched without derivatives.

The expectation can be requested with another `approx` in
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
and this check does not reject it: `"bartlett"` sums or integrates over
the support and gives a genuine expected information, at a much higher
cost per evaluation.

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
