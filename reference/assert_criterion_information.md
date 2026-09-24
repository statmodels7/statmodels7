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
answers `TRUE` and nothing is asked here. Where it does not – the
Poisson-inverse gaussian, the skew normal, the skew t and the
pseudo-Huber among the shipped families – the family falls back on the
strategy `approx` names, and the default, `"opg"`, returns \\s_i
s_i^\top\\ with \\s_i\\ the score at \\y_i\\. That is an unbiased
estimate of \\\mathcal{I}(\theta_i)\\ from one observation, not the
information itself: it depends on the response, as the observed
curvature does, and it is neither of them. The criterion built on it is
neither the Laplace approximation nor its Fisher variant, and it has no
exact outer gradient, so its search is derivative-free and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
answers `unknown`.

Measured at 4000 observations with a smooth on the first parameter, all
six families fit on the observed information in 3 to 5 criterion
evaluations with a Newton search and an analytic certificate, where the
outer-product route took 12 to 22 evaluations of a simplex. Its outer
Hessian agrees with a difference of the exact gradient at a polished
mode to between 3.8e-08 and 7.3e-08 on five of them, converging as
\\h^2\\; on the skew t it stops near 1e-06, its derivatives in \\\nu\\
being single stencils.

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
try(statmodels7:::assert_criterion_information(
  distributions7::pig1_distrib(), reml("expected"), "opg"))
#> Error : reml(hessian = "expected") reads the expected information, which 'poisson-inverse gaussian' does not write out: with approx = "opg" it would read the outer product of the scores at the data, which depends on the response and is not an expectation, so the criterion would be neither the Laplace approximation nor its Fisher variant and would have no exact outer derivatives. Use reml(hessian = "observed"), whose outer gradient and Hessian are exact on this family, or ask for the expectation by name with iwls(approx = "bartlett"), which sums or integrates over the support at a cost of minutes per evaluation.
```
