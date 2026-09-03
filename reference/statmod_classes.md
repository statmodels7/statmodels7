# The Covariance Classes of a Specification

The groups of terms that share a covariance block: those carrying the
same label and the same grouping variable, collected across the
equations.

## Usage

``` r
statmod_classes(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists: one element per
  distribution parameter, in the family's order.

## Value

A list, one element per class, each with `tag`, `group` (the deparsed
expression), `levels`, `m`, `dim` (the class total), `penalty` and
`pieces` – one per member, with its parameter, its term's name and its
within-group column count, in the order the equations were walked. Empty
where no term carries a label.

## What identifies a class

A label and a grouping, both. The label comes from
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
and the grouping from
[`modelterms7::term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.html),
which returns the expression, the levels and the within-group column
count. Two terms belong together only if the levels agree as well as the
expression: `droplevels(id)` and `id` are different expressions for one
grouping, and `id` under two subsets is one expression for two
groupings.

Correlating effects on different groupings means nothing – they are
indexed by different things and there is no block to estimate – so a
label used on two groupings is an error naming both.

## Whose prior it is

The joint prior belongs to the class and not to any of its members. At
most one term may name a `distrib`, and its dimension must be the
class's total; where none does, the default is a centered multivariate
Gaussian on
[`parameters7::dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.html),
whose coordinates are the log standard deviations and the correlations'
angles, so a printed hyperparameter is the quantity it names. At a total
of one column there is no correlation and the default is the centered
univariate Gaussian a single unlabelled term would have built.

A class whose members carry priors of **different families** – a
Gaussian intercept correlated with a Student t slope – is a copula and
not an elliptical family; modelterms7 rejects a univariate `distrib` on
a labelled term for that reason, and this function never sees the case.

## See also

[`statmod_penalty_keys()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_keys.md),
which turns each into one penalized unit;
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
for the interleaved index it is read at.
