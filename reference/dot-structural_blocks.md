# The Model's Derivative Pieces for a Filter's Recursion

Builds the `blocks` callback
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
and
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
take, at a direction or without one.

## Usage

``` r
.structural_blocks(params, ap, Vs, H, D3, D4, n)
```

## Arguments

- params:

  The distribution's parameter names.

- ap:

  Which of them carries the filter.

- Vs:

  The static rows.

- H, D3, D4:

  The family's derivatives at the fitted predictors.

- n:

  The number of observations.

## Value

A function of the direction returning a `blocks` callback.

## Details

The pieces are built on the active set the term asks for, so a panel's
outer products are of the same size whatever the number of groups.
`dcurv` serves twice, as the derivative of the curvature along the
direction and as the factor multiplying the movement of \\V_p\\, and `N`
is where the family's fourth derivative enters: each order of
differentiating the predictor through the recursion pulls in one more
order of the family.
