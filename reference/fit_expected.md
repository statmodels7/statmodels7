# Which Information Matrix a Fit Reports

`TRUE` when
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
should invert the expected information: when the fit itself inverted it,
which only
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
does and only unless asked otherwise, AND the family writes that
information out in closed form.

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

The default follows the fit rather than choosing for itself, so a
standard error comes from the same matrix the step did, and a caller who
wants the other one asks for it.

The second condition is the part that is not about the fit. Where a
family has no closed expected information the scoring step is driven by
an approximation of it – the outer product of the observed scores, which
costs one gradient and is positive semidefinite by construction – and
that is a good matrix to take a step with, the score being exact, but
not a good one to read a standard error off. Measured on a
Poisson-inverse gaussian regression at \\n = 500\\: the outer product
gives standard errors 5.7 per cent from those of the exact expectation
where the observed Hessian gives 0.6 per cent, for coefficients agreeing
to \\10^{-6}\\. So the report falls back to the observed information
there, which every family has and which is exact.

## See also

[`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html),
the predicate this reads, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which follows it.
