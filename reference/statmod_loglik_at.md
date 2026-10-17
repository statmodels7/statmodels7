# The Weighted Log-Likelihood of a Specification at Given Coefficients

Computes \$\$\ell(\beta) = \sum_i w_i \log f(y_i; \theta_i),\$\$ the
weighted log-likelihood of the whole model at one set of coefficients.
The weights enter as supplied and are never normalized, so a weight of
two counts an observation twice.

## Usage

``` r
statmod_loglik_at(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors, one per distribution parameter.

- design:

  The design. The function rebuilds it from `spec` when it is absent,
  which is convenient at the console and wasteful in a loop.

## Value

A single number. The value is `-Inf` where the density is zero at some
observation and `NaN` where the parameters are outside the family's
support; a search reads both as an unusable point.

## Details

The log-likelihood contains no penalty;
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
gives the other half of the objective, and
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
combines them.

For a model with a term of the likelihood shape, a mixture over latent
states, the log-likelihood of the term replaces the sum of the
log-densities at the observations.

## See also

[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
for its derivative,
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
for the penalized objective that it enters.
