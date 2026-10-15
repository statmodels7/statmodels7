# The Score of the Weighted Log-Likelihood

Computes the score of the weighted log-likelihood in the coefficients,
one block per distribution parameter:

\$\$\frac{\partial \ell}{\partial \beta_k} = X_k'\\(w \odot g_k)\$\$

with \\g_k\\ the per-observation derivative of the log-density in the
link-scale predictor \\\eta_k\\, which distributions7 supplies already
chained onto that scale.

## Usage

``` r
statmod_score_at(spec, coef, design = statmod_design(spec), index = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design, refreshed at `coef` if any term needs it.

- index:

  `NULL` for every coefficient, or integer positions in the stacked
  coefficients: only those columns enter the final products, which is
  what a fit over a few coordinates of a wide model reads.

## Value

A named list of numeric vectors, one per distribution parameter in the
family's order, each as long as that equation's design is wide. With
`index`, a numeric vector over those positions, in the order given.

## Details

This is the score of the likelihood with no penalty in it. The vectors
\\g_k\\ are those of
[`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md),
which already contain the contribution of a structural term. Where a
filter drives an equation, a coefficient also reaches the filter's level
through the scores at earlier times, and the reverse recursion of
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html)
adds that part. Where the likelihood is a mixture over latent states,
the score is the posterior-weighted average of the ordinary scores. The
score in the parameters of a structural term itself is returned by
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md).

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
for the value,
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
for the curvature,
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)
for a filter's own parameters.
