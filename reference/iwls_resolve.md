# Settle the Curvature of an `iwls()` Against a Family

Turns an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
left at `hessian = "auto"` into the method a fit runs: the expected
information where the family writes it exactly, and otherwise the
observed information, with the expected one taking the step where the
observed cannot. A method whose curvature was named is returned as it
is.

## Usage

``` r
iwls_resolve(method, distrib)
```

## Arguments

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- distrib:

  The distribution the fit is of.

## Value

An
[`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
object whose `hessian` is `"expected"` or `"observed"`, with `fallback`
`TRUE` exactly where `"auto"` was settled on the observed information.

## Details

[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is the question. Where it answers `FALSE` the expected information is an
approximation, the outer product of the scores by default. Measured on
68 fits of ten such families, a smooth at n = 300 and 1000: a scoring
step on that approximation ends at a point the certificate accepts in
63; the observed information alone in 60, 7 of the others ending in an
error; and the observed information with the expected standing in where
the observed penalized information is not positive definite, or where
its step finds no acceptable point, in 68 of 68, at a median time
between 0.17 and 0.98 of the expected route's per family and a REML
criterion within \\\[-3.7 \times 10^{-5}, 1.9 \times 10^{-4}\]\\ of it.
The expected information stood in 52 times over those fits, every time
for a curvature that was not positive definite.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which reads `fallback`.
