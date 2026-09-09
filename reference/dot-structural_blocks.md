# The Model's Derivative Pieces for a Filter's Recursion

Builds the `blocks` callback
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
and
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
take, at a direction or without one.

## Usage

``` r
.structural_blocks(params, ap, Vs, H, D3, D4, n, D5 = NULL)
```

## Arguments

- params:

  The distribution's parameter names.

- ap:

  Which of them carries the filter.

- Vs:

  The static rows.

- H, D3, D4, D5:

  The family's derivatives at the fitted predictors. `D5` is needed only
  where two directions are given.

- n:

  The number of observations.

## Value

A function of the direction – `NULL`, one vector, or a list of two –
returning a `blocks` callback.

## Details

The pieces are built on the active set the term asks for, so a panel's
outer products are of the same size whatever the number of groups.
`dcurv` serves twice, as the derivative of the curvature along the
direction and as the factor multiplying the movement of \\V_p\\, and `N`
is where the family's fourth derivative enters: each order of
differentiating the predictor through the recursion pulls in one more
order of the family.

With TWO directions the callback serves
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
instead, and carries three quantities more: `Q`, the fourth derivative
with two of its indices on the filter's own equation; `P`, the FIFTH
derivative contracted against both directions, which is the only place
that order enters; and `cppp`, the scalar the level's own curvature
moves by, which no contraction recovers. `N` is then a list of two, one
per direction.
