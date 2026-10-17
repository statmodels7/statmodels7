# The Linear Predictors and the Parameters That They Give

Turns a coefficient structure into the per-observation parameters that
every distributions7 generic takes. Each equation's design,
coefficients, offset and adjustment give a linear predictor \\\eta_k\\,
and the parameter's own link inverts it to \\\theta_k\\.

## Usage

``` r
statmod_eta(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  read for the distribution, its links, the offsets and the sample size.

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it, refreshed at `coef` if any term needs it.

- coef:

  A named list of coefficient vectors, one per distribution parameter,
  each as long as its equation's design is wide.

## Value

A list with the elements

- `eta`:

  a named list of the linear predictors, one numeric vector of length
  `spec@n_obs` per distribution parameter in the family's order.

- `theta`:

  the parameters that the predictors imply, in the same shape.

- `filters`:

  the runs of the filter-shaped structural term, an empty list when the
  model has none.

- `regimes`:

  the latent-state structure of a term of the likelihood shape, an empty
  list when the model has none.

- `eta_static`:

  present only when the model carries a structural term: the predictors
  before the term's level is added.

## Details

\$\$\eta_k = X_k \beta_k + o_k + a_k, \qquad \theta_k =
g_k^{-1}(\eta_k)\$\$

with \\o_k\\ the equation's offset and \\a_k\\ the adjustment that a
refreshable term contributes, zero for an ordinary design.

The inverse link is evaluated by
[`linkfunctions7::linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.html).
The links onto a bounded interval keep the result within the bounds of
the parameter, at or next to a bound for a predictor far out on the link
scale, so that a density can still be evaluated. The `inverse` link is
an exception: a negative predictor gives a negative parameter.

A structural term of the filter shape adds its level to \\\eta_k\\
before the link is inverted. A term whose contribution is a likelihood
mixed over latent states reports the posterior-weighted predictor.

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
and
[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md),
which read this,
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the refresh.
