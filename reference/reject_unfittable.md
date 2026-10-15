# Reject a Term the Fitting Scheme Does Not Cover

Signals an error naming any term that the fitting scheme of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does not cover, and any formula that carries more than one structural
term.

## Usage

``` r
reject_unfittable(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists, one per distribution
  parameter.

## Value

`NULL`, invisibly; called for the error.

## Details

Terms are examined by what they implement and not by a list of class
names, so a term written later is covered without an edit here.

A **structural** term rewrites the likelihood or the filtered predictor
instead of contributing a fixed design block, so it has no design block
and answers neither `term_matrix()` nor `term_npar()`. A term of the
filter shape, such as `gas()`, contributes a predictor through
`term_filter()`. A term of the likelihood shape, such as `regime()`,
contributes a likelihood through `term_loglik()`.
[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
routes both. What this function rejects is a structural term that
implements neither `term_filter()` nor `term_loglik()`, and a structural
term of the likelihood shape that carries a covariance label
([`unfittable_reason()`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md)).
[`reject_incompatible()`](https://statmodels7.github.io/statmodels7/reference/reject_incompatible.md)
is called before these reasons are raised, so a formula with more than
one structural term is rejected with that function's message.

A term whose block depends on its own coefficients is fitted by
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
and is not rejected here.

## See also

[`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
