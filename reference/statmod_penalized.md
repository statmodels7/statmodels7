# Every Penalized Unit of a Specification

Returns one entry per penalty in the model, regardless of the term to
which the penalty belongs and of whether that term has more than one
penalty.

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

The key is the term's name in the formula. When the term carries more
than one penalty and the entry is named, the entry's own name is
appended after `::`. Two `ridge()` terms are two terms with two keys and
two hyperparameters. A term with one penalty has only its name as key,
so a hyperparameter is read by the term's name.

`statmod_penalized_units()` builds the list. `statmod_penalized()`
returns the list from the design's `eta_memo` environment when the
design has one and the terms and the family's parameters are the ones
for which the design was built. An objective evaluated many times over
one design reads the list from there.

## See also

[`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md),
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
