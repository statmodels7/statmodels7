# Let Every Term Draw the Coefficients Only It Can

Walks the terms of every equation and hands each the slice of the drawn
coefficients its block owns, taking back whatever it chose to replace.

## Usage

``` r
rstatmod_term_draw(spec, design, coef, params, fixed_eq, sd)
```

## Arguments

- spec:

  The specification, read for its terms.

- design:

  The design, read for each term's columns.

- coef:

  The coefficients drawn so far, a list by parameter.

- params:

  The distribution parameter names.

- fixed_eq:

  The parameters `par` fixes whole.

- sd:

  The width of the draws.

## Value

A list with `coef`, the coefficients with each term's own replaced, and
`owned`, a list by parameter of named numeric vectors carrying one width
per replaced position, the names being the positions.

## Details

Almost every term replaces nothing: a slope is measured in the
response's units against a covariate's and a normal of the caller's
width is as good a truth as any. A break-point is the exception, being a
position on the covariate's own axis, and
[`modelterms7::term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.html)
is where a term says so.

It runs last, after the plain draw and after the priors, because a
coordinate a term owns is one no other rule describes correctly. Which
coordinates those were comes back so that a prior covering them is
reported at the width the term used rather than at the one it drew.

An equation `par` fixes whole is skipped: its coefficients are the
caller's and are written over everything afterwards anyway.

The walk is over the design's blocks and does not descend into a term's
subformulas, so a break-point developing another term's own parameter –
`nl(a ~ 0 + seg(x))` – is drawn by the plain rule as it was before. The
sub-terms a break-point's own development carries are `linpar()` and
`random()`, which have nothing to say here, so the shapes that reach one
are covered by the block walk.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md),
[`modelterms7::term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.html)
