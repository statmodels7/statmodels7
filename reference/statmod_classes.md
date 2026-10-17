# The Covariance Classes of a Specification

Returns the groups of terms that share a covariance block: those
carrying the same label and the same grouping variable, collected across
the equations.

## Usage

``` r
statmod_classes(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists: one element per
  distribution parameter, in the family's order.

## Value

A list named by the class key (`<label> | <grouping>`), one element per
class, each with `tag`, `group` (the deparsed expression), `levels`,
`m`, `key`, `dim` (the class total), `space` (see
[`class_space()`](https://statmodels7.github.io/statmodels7/reference/class_space.md)),
`structural`, `sterm` (the name of the structural term, present only
where the class touches one), `penalty` and `pieces`. The pieces are one
per member, as
[`label_pieces()`](https://statmodels7.github.io/statmodels7/reference/label_pieces.md)
returns them, in the order in which the equations were walked. The list
is empty if the terms carry no label.

## What identifies a class

A class is identified by both a label and a grouping. The label comes
from
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
and the grouping from
[`modelterms7::term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.html),
which returns the expression, the levels and the within-group column
count. Two terms belong together only if the levels agree as well as the
expression: `droplevels(id)` and `id` are different expressions for one
grouping, and `id` under two subsets is one expression for two
groupings.

Correlating effects on different groupings means nothing, since they are
indexed by different things and there is no block to estimate, so a
label used on two groupings is an error naming both.

## Whose prior it is

The joint prior belongs to the class and not to any of its members. At
most one term may name a `distrib`, and its dimension must be the
class's total; if no term names one, the default is a centered
multivariate Gaussian on
[`parameters7::dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.html),
whose coordinates are the log standard deviations and the correlations'
angles, so a printed hyperparameter is the quantity that it names. At a
total of one column there is no correlation and the default is the
centered univariate Gaussian that a single unlabeled term would have
built.

A class whose members carry priors of different families (a Gaussian
intercept correlated with a Student t slope) is a copula and not an
elliptical family. modelterms7 rejects a univariate `distrib` on a
labeled term for that reason, so the case does not reach this function.

## See also

[`statmod_penalty_keys()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_keys.md),
which turns each into one penalized unit;
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
for the interleaved index at which it is read.
