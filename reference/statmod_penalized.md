# Every Penalized Unit of a Specification

One entry per penalty in the model, whatever term it belongs to and
whether or not that term has more than one.

## Usage

``` r
statmod_penalized(spec, design)

statmod_penalized_units(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A list of entries, each with `param`, `term` (the name in the formula),
`key`, `within`, `penalty`, `fixed`, `n_values`, `values`, `min_ratio`,
`search`, `ids`, `structural`, `cols` and `index`. For an ordinary unit
`cols` are positions within the parameter's coefficients and `index`
positions in the stacked vector. For a structural unit `cols` are
positions among the term's own parameters and `index` is `NULL`. The
entry of a covariance class also carries `class`, `params`, `pieces`,
`mixed` and, for a mixed class, `joint` and `beta_index`.

## Details

The enumeration runs over the distribution parameters and over the terms
of each, and reads the penalties of each term. A term may carry several
penalties, over different subsets of its own parameters, as a panel
model with a population value and a shrunk deviation per group needs. A
covariance class (see
[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md))
is one unit over the stacked columns of its members. The callers that
need the penalties of a model read this one enumeration.

**The key** is the term's name in the formula, with the entry's own name
appended after `::` when the term carries more than one penalty and the
entry is named. Two `ridge()` terms are two terms with two keys and two
hyperparameters. A term with one penalty has its name alone as key, so a
hyperparameter is read by the term's name.

`statmod_penalized_units()` builds the list. `statmod_penalized()`
returns it from the design's `eta_memo` environment where the design has
one and the terms and the family's parameters are the ones it was built
for, which is what an objective evaluated many times over one design
reads.

## See also

[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md),
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
