# The Linear Predictors and the Parameters They Give

Turns a coefficient structure into the per-observation parameters every
distributions7 generic takes. Each equation's design, coefficients,
offset and adjustment give a linear predictor \\\eta_k\\, and the
parameter's own link inverts it to \\\theta_k\\.

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

A list of two named lists, each with one entry per distribution
parameter in the family's order:

- `eta`:

  the linear predictors, numeric vectors of length `spec@n_obs`.

- `theta`:

  the parameters they imply, the same shape, each strictly inside its
  parameter's bounds.

## Details

\$\$\eta_k = X_k \beta_k + o_k + a_k, \qquad \theta_k =
g_k^{-1}(\eta_k)\$\$

with \\o_k\\ the equation's offset and \\a_k\\ the adjustment a
refreshable term contributes, zero for an ordinary design.

The inverse link clamps \\\theta_k\\ strictly inside its own open
bounds, so a predictor far out on the link scale gives a parameter a
density can be evaluated at. A structural term of the filter shape adds
its level to \\\eta_k\\ before the link is inverted.

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
and
[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md),
which read this,
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the refresh.
