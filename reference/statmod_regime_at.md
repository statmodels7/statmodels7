# Run the Structural Terms at the Current Parameters

The level each filter adds to the equation it sits in, together with the
derivative of that level in the term's own parameters.

The smoothed state probabilities of each latent Markov term, the levels
its regimes shift the predictor by, and the per-state predictors those
give.

## Usage

``` r
statmod_regime_at(spec, design, eta_static, theta_static)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- eta_static:

  The static predictors.

- theta_static:

  The parameters they give.

## Value

A named list, one entry per structural term, each with `param`, `level`,
`jacobian` and the callbacks used.

A list, one entry per regime term.

## Details

The result is memoized on the coefficients and the term's parameters,
since the objective, its gradient and its curvature are asked for at the
same point in turn and a filter is the expensive part of each.

A term of this shape does not report a predictor: its contribution is a
likelihood mixed over states, so the model's log-likelihood comes from
[`modelterms7::term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.html)
and never from the density at a single predictor.

Everything else follows from Fisher's identity: the derivative of that
mixed likelihood in any predictor the model carries is the
posterior-weighted derivative of the ordinary one. A caller therefore
differentiates its own log-density once per state, vectorized, and
weights.

The predictor reported for the equation the term sits in is the
posterior-weighted one, and that is what a fitted value means here.

## See also

[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)

[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
