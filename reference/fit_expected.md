# Which Information Matrix a Fit Reports

Returns `TRUE` when
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
inverts the expected information: when the fit itself inverted it, which
only
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
does and only with `hessian = "expected"`, and when the family writes
that information out in closed form.

## Usage

``` r
fit_expected(object)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single logical.

## Details

The default follows the fit, so a standard error comes from the same
matrix that the step used, and a caller obtains the other matrix by
setting `expected` explicitly.

The second condition does not concern the fit. Where a family has no
closed expected information, the scoring step is driven by an
approximation of it, the outer product of the observed scores, which
costs one gradient and is positive semidefinite by construction. That
matrix is adequate for taking a step, because the score is exact, but it
is a poor basis for a standard error, which it can place far from the
one given by the exact expectation. The report therefore falls back to
the observed information there, which every family has and which is
exact.

## See also

[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html),
the predicate that this function reads, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which follows it.
