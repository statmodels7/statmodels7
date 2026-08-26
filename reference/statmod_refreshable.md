# Which Terms Recompute Their Own Block

Locates every term of the specification whose design block is a function
of its own coefficients, and reports where each one sits. These are the
terms whose block has to be rebuilt whenever the coefficients move:
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
and
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html).
Every other term's block is fixed once at build time.

## Usage

``` r
statmod_refreshable(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  whose `terms` are walked equation by equation.

## Value

A list with one element per refreshable term, in the order the design
holds them: equations in the family's parameter order, and terms within
an equation in the order they were written. Each element is a list of
two:

- `param`:

  the distribution parameter whose equation the term sits in, a string.

- `term`:

  the term's index within that equation, an integer.

An empty list when no term refreshes, which is the common case.

## Details

Membership is decided by asking whether the term registers a
`term_refresh()` method of its own, which is
[`refreshes_own_block()`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md).
The base method on `model_term` is the identity, so a term written later
is covered without an edit here.

This is the list every other function in the file walks. Its emptiness
is what makes a model of ordinary terms pay nothing for the refresh
machinery.

## See also

[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the rebuild,
[`refreshes_own_block()`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md)
for the predicate,
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
for advancing the state.
