# The Objective, Its Gradient and Its Hessian, Stacked

\\F(\beta) = -\ell(\beta) + \sum_t \rho_t\\, unaveraged, over the
coefficients of every parameter stacked into one vector.

## Usage

``` r
statmod_objective(
  spec,
  hyper,
  design = statmod_design(spec),
  expected = TRUE,
  approx = "opg"
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- hyper:

  The hyperparameters, per penalized term, held for the life of the
  objective.

- design:

  The design. Rebuilt from `spec` when absent.

- expected:

  `TRUE` for the expected information in `he`, `FALSE` for the observed
  one.

- approx:

  How the expected information is approximated for a family with no
  closed form.

## Value

A list of five functions:

- `fn(b)`:

  the objective at the stacked coefficients `b`, a single number.

- `gr(b)`:

  its gradient, a numeric vector as long as `b`.

- `he(b)`:

  its Hessian, a `p x p` matrix.

- `split(b)`:

  the stacked vector as a named list, one entry per distribution
  parameter.

- `stack(l)`:

  the inverse of `split`.

## The objective is not divided by the sample size

A penalty is a negative log-prior at full size, and a posterior adds a
log-likelihood and a log-prior at full size. Averaging the likelihood
alone would make a hyperparameter mean something that depends on \\n\\.

What is scaled instead is the stopping rule, in the one place it is
read. See
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
and
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md).

## The five closures share one design

All five close over `design`, `hyper` and the two information settings,
so a caller who changes a hyperparameter builds a new objective. That is
what the outer search does at every point it visits.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
and
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which consume this,
[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
and
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
for its two halves.
