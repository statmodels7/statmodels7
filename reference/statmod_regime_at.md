# Run the Terms of the Likelihood Shape at the Current Parameters

The smoothed posterior over the components of each structural term of
the likelihood shape (the states of a regime term, the intervals or
nodes of a marginal break-point term), the shifts that the components
add to the predictor, and the posterior-weighted predictor that these
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

A named list, one entry per structural term of the likelihood shape (a
regime term or a marginal break-point term), each with `param`, `term`,
`tm`, `psi`, `gamma` (the posterior), `mu` (the shifts), `cb` (the
callbacks), `eta_static` and `eta_bar` (the posterior-weighted
predictor).

## Details

A term of this shape does not report a predictor: its contribution is a
likelihood mixed over components, so the log-likelihood of the model
comes from
[`modelterms7::term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.html)
and not from the density at a single predictor.

The remaining quantities follow from Fisher's identity: the derivative
of that mixed likelihood in any predictor of the model is the
posterior-weighted derivative of the ordinary one. A caller therefore
differentiates its own log-density once per component, vectorized, and
weights.

The predictor reported for the equation the term sits in is the
posterior-weighted one, and that is what a fitted value means here.

## See also

[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
