# The Score of the Weighted Log-Likelihood

Computes the score of the weighted log-likelihood in the coefficients,
one block per distribution parameter:

\$\$\frac{\partial \ell}{\partial \beta_k} = X_k'\\(w \odot g_k)\$\$

with \\g_k\\ the per-observation derivative of the log-density in the
link-scale predictor \\\eta_k\\, which distributions7 supplies already
chained onto that scale.

## Usage

``` r
statmod_score_at(spec, coef, design = statmod_design(spec))
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

A named list of numeric vectors, one per distribution parameter in the
family's order, each as long as that equation's design is wide.

## Details

This is the score of the likelihood alone, with no penalty in it, and it
is the gradient of a model whose design is fixed.

A model carrying a structural term needs a correction that is not here:
the term's level is driven by scores read at earlier predictors, so a
coefficient reaches it through the recursion.
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)
and
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html)
supply that.

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
for the value,
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
for the curvature,
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)
for a filter's own parameters.
