# Every Labelled Effect a Term Carries, Its Sub-Terms Included

One piece per labelled random effect reachable from a term: the term
itself where it carries a label, and the sub-terms developing its own
parameters otherwise, walked to any depth.

## Usage

``` r
label_pieces(
  term,
  param,
  nm,
  within = NULL,
  path = character(0),
  structural = NULL
)
```

## Arguments

- term:

  One built term.

- param:

  The distribution parameter its equation belongs to.

- nm:

  The equation-level term's name, which is what the design is keyed by.

- within:

  The piece's columns in that term's block, or `NULL` at the top level,
  where the piece is the whole of it.

- path:

  Which of the parent's parameters were developed to reach the piece,
  outermost first, empty at the top level. It says what the effect is an
  effect ON: a labelled random intercept written inside
  `seg(x, psi ~ ...)` is an effect on the break-point and not on the
  mean, and a report naming only the equation would be read as the
  second.

- structural:

  Whether `within` indexes a structural term's own parameters rather
  than columns of the design. `NULL` at the top level, where it is read
  from the term; a recursive call passes what it was told.

## Value

A list of pieces, each with `param`, `term`, `within`, `path`, `dim`,
`tag`, `group` (as
[`modelterms7::term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.html)
returns it), `structural` and `distrib`. Empty where nothing under the
term is labelled.

## Why a sub-term is reachable at all

`seg(x, psi ~ random(~ 1 | u | id))` develops a break-point over a
labelled random effect. The labelled term is not one of the equation's
terms – the equation carries one `SegTerm`, whose own label is absent –
and its coefficients are columns of that term's block. Measured, they
are exactly the ones
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html)
reports as that component's `sub_index`, so a piece records them as
`within`, positions in the parent's block, and the design turns them
into positions in the stacked vector the same way it does for any other
term.

That is the whole of what a subformula costs here, and it is why the
case is covered: a labelled effect written in a subformula of an
**additive** term lives in the same vector as one written in an
equation.

A **structural** parent addresses its coefficients differently: they are
the term's own parameters, which contribute no design column and live in
the design's structural state. The walk is the same and the positions it
composes are the same numbers – measured, a labelled effect inside
`gas(alpha1 ~ 1 + random(~ 1 | u | g))` comes out at 3 to 12, exactly
the `cols` the unlabelled sub-term's own penalty is read at – so what a
piece records is which vector they index rather than a different
arithmetic.

## Depth

A sub-term is an ordinary term and may develop parameters of its own, so
the walk recurses and composes `within` on the way down: a depth-two
sub-term's columns are positions in its parent's block, which are
themselves positions in the equation-level term's.

A labelled term's own sub-terms are not walked. Its columns are the
class's already, and anything inside it belongs to that block rather
than to another.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md),
its caller;
[`class_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_pieces.md)
for the mapping of `within` onto the stacked vector.
