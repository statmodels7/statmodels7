# The Covariance Labels a Term Carries, Its Sub-Terms Included

Every covariance label reachable from a term: its own, and those of the
sub-terms developing its own parameters, walked to any depth.

## Usage

``` r
term_tags_deep(term, what = "it")
```

## Arguments

- term:

  One built term.

- what:

  How to name this term in the result's names. Defaults to the term
  itself, and a recursive call passes the path it came by.

## Value

A named character vector, empty where no label is reachable.

## Details

A label written inside a subformula is not visible where the equation's
terms are enumerated: `seg(x, psi ~ random(~ 1 | u | id))` is one
`SegTerm`, whose own
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
is `NA`, and the labelled term is a sub-term of its break-point. The
sub-terms are reached through
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html),
whose `subs` field is the list of them, and the walk recurses because a
sub-term is an ordinary term and may develop parameters of its own.

The names of the result say where each label was found, so a message can
name the sub-term rather than the equation's term alone.

## See also

[`unfittable_reason()`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md),
its caller;
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
for one term's answer.
