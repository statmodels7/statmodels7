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
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A named list, one character vector per structural term.

## Details

A score-driven level and a regime's first level both add a constant to
their equation's predictor. With an intercept there too the two are
exactly confounded: shifting the intercept by \\c\\ and the level by
\\-c(1 - \sum_j \xi_j)\\ leaves every predictor unchanged, since the
recursion is affine in the level given the score path, the score path
depends on the predictor alone, and the starting level moves by the same
\\c\\. The likelihood is flat along that direction, and a fit can reach
the ridge without failing, because the score is small wherever the
surface is flat.

**The linear intercept wins.** Where both are present the term's level
is held at zero and the coefficient carries it, which makes
`y ~ x + gas(...)` an ordinary formula. The model loses nothing, since a
constant cannot express the dynamics or the difference between one
regime and another, and the parameters that express them remain free.

The test is made on the **span** of the equation's design and not on a
column named `"(Intercept)"`: a factor coded without one, or any set of
columns summing to a constant, spans the constant just as well.
[`modelterms7::term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.html)
gives which parameter of the term is the level.

**For a developed level** the same test applies to a subspace. With
`omega ~ Z gamma` the confounding is no longer with one constant but
with whatever `span(Z)` shares with the span of the equation's design.
The constant coordinates are held as above. For the rest, an unpenalized
coordinate whose column lies in the equation's span is flagged with a
warning and is not held: holding it would change the model where the
confounding is not exact (a time-varying shared column is exactly flat
only when its lags stay in the development's span, which depends on
\\q\\ and on the column). A penalized coordinate is identified by its
penalty. Where the direction really is flat, the variance matrix shows
it.

## See also

[`reject_incompatible()`](https://statmodels7.github.io/statmodels7/reference/reject_incompatible.md)
