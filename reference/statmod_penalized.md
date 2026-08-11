# Every Penalized Unit of a Specification

One entry per penalty in the model, whatever term it belongs to and
whether or not that term has more than one.

## Usage

``` r
statmod_penalized(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A list of entries, each with `param`, `term` (the name in the formula),
`key`, `cols` (positions within the parameter's coefficients), `index`
(positions in the stacked vector) and `penalty`.

## Details

Twelve places used to run the same loop – over the distribution
parameters, over each one's terms, asking each term for its penalty –
and each of them assumed a term carries at most one. A term may carry
several, over different subsets of its own parameters, which is what a
panel model with a population value and a shrunk deviation per group
needs. Enumerating once is both the generalization and the removal of
eleven copies.

**The key** is the term's name in the formula, and the entry's own name
appended after `::` when the term carries more than one. Two `ridge()`
terms are two terms with two keys and two hyperparameters, which they
already were; a term with one penalty over the whole of itself keys
exactly as before, so nothing that reads a hyperparameter by term name
changes.

## See also

[`statmod_blocks`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md),
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
