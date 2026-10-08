# Reject Two Terms of One Equation That Name Their Coefficients Alike

The coefficients of an equation are addressed by name in
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) and the summary, so
two terms of one equation must not give their coefficients the same
names. A term prefixes its names with its label, which for `s()` and
`te()` is built from the covariates: two smooths of the same covariate
with different constructions, such as
`s(x, bspline_smooth(k = 10)) + s(x, pspline_smooth())`, would share it.
Such a formula is rejected here, naming the terms and the argument
`label` that separates them.

## Usage

``` r
reject_duplicate_names(terms)
```

## Arguments

- terms:

  The built terms, a list with one element per distribution parameter,
  each a named list of terms as
  [`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
  builds them.

## Value

`NULL`, invisibly. Called for the error.
