# The Objective, Its Gradient and Its Hessian, Stacked

Builds the objective \\F(\beta) = -\ell(\beta) + \sum_t \rho_t\\, which
is not averaged over the observations, over the coefficients of every
parameter stacked into one vector.

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

  The design. The function rebuilds it from `spec` when it is absent.

- expected:

  `TRUE` for the expected information in `he`, `FALSE` for the observed
  one.

- approx:

  How the expected information is approximated for a family with no
  closed form.

## Value

A list with the elements

- `npar`:

  the number of coefficients of each equation, an integer vector.

- `fn(b)`:

  the objective at the stacked coefficients `b`, a single number.

- `gr(b)`:

  its gradient, a numeric vector as long as `b`.

- `gr_sub(b, idx)`:

  the gradient only at the positions `idx`.

- `he(b)`:

  its Hessian, the information matrix (expected or observed according to
  `expected`) plus the Hessian of the penalties, a `p x p` matrix.

- `split(b)`:

  the stacked vector as a named list, one entry per distribution
  parameter.

- `stack(l)`:

  the inverse of `split`.

## The objective is not divided by the sample size

A penalty is a negative log-prior at full size, and a posterior adds a
log-likelihood and a log-prior at full size. Averaging only the
likelihood would make a hyperparameter mean something that depends on
\\n\\.

Instead, the stopping rule is scaled, in the one place where it is read.
See
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
and
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md).

## The closures share one design

The functions `fn`, `gr`, `gr_sub` and `he` close over `design` and
`hyper`, and `he` also closes over the two information settings, so that
a change of a hyperparameter requires a new objective. The outer search
builds a new objective at every point that it visits.

## See also

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
and
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which consume this,
[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
and
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
for its two halves.
