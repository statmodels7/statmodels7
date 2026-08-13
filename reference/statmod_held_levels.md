# Which Structural Levels a Linear Intercept Already Carries

The parameters of the structural terms that must be held rather than
estimated, because the equation they sit in already spans the constant
they would shift it by.

## Usage

``` r
statmod_held_levels(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A named list, one character vector per structural term.

## Details

A score-driven level and a regime's first level both add a constant to
their equation's predictor. With an intercept there too the two are
EXACTLY confounded: shifting the intercept by \\c\\ and the level by
\\-c(1 - \sum_j b_j)\\ leaves every predictor unchanged, since the
recursion is affine in the level given the score path, the score path
depends on the predictor alone, and the starting level moves by the same
\\c\\. The likelihood is flat along that direction, and a fit reaches
the ridge without failing – the score is small because the surface is
flat, not because it is a maximum.

**The linear intercept wins.** Where both are present the term's level
is held at zero and the coefficient carries it, which is what makes
`y ~ x + gas(...)` an ordinary thing to write. Nothing about the model
is lost: what a constant cannot express is the dynamics, or the
difference between one regime and another, and those are the parameters
that remain free.

The question is asked of the **span** of the equation's design and not
of a column named `"(Intercept)"`: a factor coded without one, or any
set of columns summing to a constant, spans it just as well. Which
parameter is the level is the term's own answer, through
[`term_level_param`](https://statmodels7.github.io/modelterms7/reference/term_level_param.html).

**A developed level asks the same question of a subspace.** With
`omega ~ Z gamma` the confounding is no longer with one constant but
with whatever `span(Z)` shares with the span of the equation's design.
The constant coordinates are the term's own answer, held as above. For
the rest, an unpenalized coordinate whose column lies in the equation's
span is flagged with a warning rather than held: holding it would change
the model where the confounding is not exact (a time-varying shared
column is exactly flat only when its lags stay in the development's
span, which depends on \\q\\ and on the column), while a penalized
coordinate is identified by its penalty, exactly as a deviation is.
Where the direction really is flat, the variance matrix names it.

## See also

[`reject_incompatible`](https://statmodels7.github.io/statmodels7/reference/reject_incompatible.md)
