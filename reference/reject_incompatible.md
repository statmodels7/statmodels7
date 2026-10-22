# Combinations of Terms That Are Not a Model

Rejects a formula carrying more than one structural term, in whatever
equations they sit.

## Usage

``` r
reject_incompatible(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists.

## Value

`NULL`, invisibly; called for the error.

## Details

A filter is driven by the score of the log-likelihood at the predictor
that it has just produced, so two filters in one equation are two levels
adding up, with nothing to tell one from the other, and are rejected.
Score-driven terms in different equations are accepted: each is driven
by a score that depends on the others' state, so together they are one
recursion, which
[`join_gas_terms()`](https://statmodels7.github.io/statmodels7/reference/join_gas_terms.md)
assembles. Any other structural term beside a second one is rejected. A
term whose contribution is a likelihood mixed over latent states does
not report a predictor at all, so it cannot be combined with anything
that expects one.

The count is over the whole formula and not per equation, because a
structural term is admissible only where everything else in the model is
a predictor that can drive it.

## See also

[`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md),
[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
