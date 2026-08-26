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

  The design. Rebuilt from `spec` when absent, which is convenient at
  the console and wasteful in a loop.

## Value

A single number. `-Inf` where the density is zero at some observation,
and `NaN` where the parameters are outside the family's support, both of
which a search reads as an unusable point.

## Details

No penalty enters. This is the likelihood alone;
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
gives the other half of the objective, and
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
combines them.

## See also

[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
for its derivative,
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
for the penalized objective it enters.
