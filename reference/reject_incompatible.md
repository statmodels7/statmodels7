# Combinations of Terms That Are Not a Model

Rejects a formula carrying more than one structural term, whatever
equations they sit in.

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

Two of them are not a model that the layer could fit and then report. A
filter is driven by the score of the log-likelihood at the predictor it
has just produced, so two filters in one equation are two levels adding
up, with nothing to tell one from the other; in two equations, each is
driven by a score that depends on the other's state, and the pair is one
recursion written as two, which neither term implements. A term whose
contribution is a likelihood mixed over latent states does not report a
predictor at all, so it cannot be combined with anything that expects
one.

The count is over the whole formula rather than per equation because
that is the honest boundary: what makes one admissible is that
everything else in the model is a predictor it can be driven by.

## See also

[`reject_unfittable`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md),
[`statmod_structural`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
