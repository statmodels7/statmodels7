# The Score of the Weighted Log-Likelihood, Observation by Observation

Computes, for every observation, the derivative of its weighted
log-likelihood in the static predictor of each equation.
[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
projects these vectors onto the columns of the design, and
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
builds from them the working response of a block that is fitted by
coordinate descent.

## Usage

``` r
statmod_score_obs(spec, coef, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design, refreshed at `coef` if any term needs it.

## Value

A named list of numeric vectors of length `spec@n_obs`, one per
distribution parameter in the family's order. Each entry already carries
the observation weights.

## Details

The static predictor of an equation is the sum of its columns times
their coefficients, with the offset and the adjustment of a block that
moves with its coefficients. In an ordinary model it is the whole
predictor, and the derivative is the family's score on the link scale
times the observation weight.

A structural term adds to the static predictor a part of its own, and
the derivative then carries more than the family's score at the whole
predictor. Where a filter drives an equation, a coefficient also reaches
the filter's level through the scores at earlier times, and the reverse
recursion of
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html)
adds that part. Where the likelihood is a mixture over latent states,
Fisher's identity gives the derivative as the posterior-weighted average
of the ordinary score over the states.

## See also

[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md),
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
