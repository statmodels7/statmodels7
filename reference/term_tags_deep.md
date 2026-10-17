# The Covariance Labels That a Term Carries, Including Its Sub-Terms

Returns every covariance label reachable from a term: its own, and those
of the sub-terms that develop its own parameters, found by walking to
any depth.

## Usage

``` r
term_tags_deep(term, what = "it")
```

## Arguments

- term:

  One built term.

- what:

  How to name this term in the result's names. Defaults to the term
  itself, and a recursive call passes the path by which the term was
  reached.

## Value

A named character vector, empty when the term has no reachable label.

## Details

A label written inside a subformula is not visible where the equation's
terms are enumerated: `seg(x, psi ~ random(~ 1 | u | id))` is one
`SegTerm`, whose own
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
is `NA`, and the labeled term is a sub-term of its break-point. The
sub-terms are reached through
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html),
whose `subs` field is the list of them, and the walk recurses because a
sub-term is an ordinary term and may develop parameters of its own.

The names of the result record where each label was found, so a message
can name the sub-term rather than only the equation's term.

## See also

[`unfittable_reason()`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md),
its caller;
[`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html)
for one term's answer.
