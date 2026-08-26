# Reject a Term the Fitting Scheme Does Not Cover

Signals an error naming any term whose block is not a fixed design,
which is what the alternation of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
assembles.

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

One shape is outside that assembly, and it is read off the term rather
than from a list of class names, so a term written later is covered
without an edit here.

A **structural** term rewrites the likelihood instead of contributing a
predictor, so it has no design block at all and answers neither
`term_matrix()` nor `term_npar()`. Reaching it through the design
produced an error naming one of those generics, which says nothing about
the cause.
[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
routes those, and what remains here is the term class that is structural
and implements neither shape of the contract.

A term whose block depends on its own coefficients was rejected here too
until the alternation learned to refresh one; it is fitted now, by
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md).

Every equation is examined before the error is raised, so a model
carrying one such term in the mean and another in the scale reports both
rather than the first.

## See also

[`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
