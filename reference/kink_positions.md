# Where the Objective Has a Kink, in the Stacked Coefficients

Asks every term that recomputes its own block where, among its own
coefficients, the objective is not differentiable at the coefficients
given, through
[`modelterms7::term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.html),
and carries the answer into the positions of the stacked vector.
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
holds those positions where a line search rejects every step length.

## Usage

``` r
kink_positions(spec, design, coef, split)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, a named list with one numeric vector per parameter.

- split:

  The stacked positions of each parameter's coefficients, a named list
  as the objective's `split` returns it for
  [`seq_along()`](https://rdrr.io/r/base/seq.html) of the stacked
  vector.

## Value

An integer vector of stacked positions, possibly empty.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).
