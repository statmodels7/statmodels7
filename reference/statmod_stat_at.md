# One Likelihood Statistic for One Coefficient

The Wald, likelihood-ratio, Rao score or Terrell gradient statistic for
the hypothesis that one coefficient equals `value`, with its p-value.

## Usage

``` r
statmod_stat_at(
  fit,
  param,
  coefname,
  value = 0,
  test = "wald",
  type = "bayesian",
  mode_error = FALSE,
  objective0 = NULL
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, coefname:

  Which coefficient, as in
  [`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md).

- value:

  The value under the null, a single number. `0` by default.

- test:

  One of `"wald"`, `"lr"`, `"score"`, `"gradient"`.

- type:

  Which variance the Wald statistic reads, as
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  takes it.

- mode_error:

  Whether to read how far above its mode the restricted fit stopped
  ([`restricted_mode_error()`](https://statmodels7.github.io/statmodels7/reference/restricted_mode_error.md)).
  `FALSE` by default, because it costs one Hessian where the likelihood
  ratio and the gradient statistic need none.
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
  which calls this once per row, and
  [`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
  which calls it once per evaluation of the statistic, do not read it.

- objective0:

  The penalized objective at the unrestricted joint mode, as
  [`joint_objective()`](https://statmodels7.github.io/statmodels7/reference/joint_objective.md)
  returns it, for a caller that reads the likelihood ratio many times on
  one fit. `NULL` computes it here.

## Value

A list with `test`, `statistic`, `df`, `p.value`, `converged` – `NA` for
the Wald statistic and the restricted fit's flag for the other three –
and `mode_error`, `NA` unless it is requested.

## Details

All four read the penalized objective, which is the function the fit
maximizes. Writing \\\ell_p(\beta) = \ell(\beta) - \rho(\beta;\theta)\\
for it, \\\hat\beta\\ for the unrestricted estimates, \\\tilde\beta\\
for those with \\\beta_j\\ held at \\b\\, \\U = \partial\ell_p\\ for its
score and \\K\\ for the penalized information, the four are \$\$W =
(\hat\beta_j - b)^2 / \widehat{\mathrm{Var}}(\hat\beta_j), \qquad LR =
2\\\ell_p(\hat\beta) - \ell_p(\tilde\beta)\\,\$\$ \$\$S =
U_j(\tilde\beta)^2 \\ \[K(\tilde\beta)^{-1}\]\_{jj}, \qquad T =
U_j(\tilde\beta)\\(\hat\beta_j - b),\$\$ each compared with a
\\\chi^2_1\\. Where the model carries no penalty \\\ell_p\\ is the
log-likelihood and these are the classical four.

The score statistic has the general form \\U'K^{-1}U\\. It reduces to
the product above because every other component of \\U\\ vanishes at the
restricted maximum, which also indicates that the restricted fit has
converged.

## The penalized objective and the likelihood

The likelihood ratio is taken on \\\ell_p\\ and not on \\\ell\\.
\\\tilde\beta\\ maximizes \\\ell_p\\ under the restriction and not
\\\ell\\, so \\\ell(\tilde\beta)\\ can exceed \\\ell(\hat\beta)\\ and
the difference of log-likelihoods can be negative. The other three
statistics read the penalized objective as well: \\U\\ and \\K\\ are its
gradient and its curvature, and the Wald variance is \\(H+S)^{-1}\\
under `type = "bayesian"`.

For a shrunk coordinate this is the reading that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
gives under `type = "bayesian"`: the penalty is a log-prior, so
\\\ell_p\\ is a log posterior at the hyperparameters in force, and the
interval by inversion is the set that the posterior does not reject,
conditional on them. The interpretation that the Wald interval on such a
row carries therefore extends to the other three statistics. On a smooth
the inverted interval and the Wald one agree to second order, and the
inversion adds the asymmetry.

## Differences between the statistics

The Wald statistic reads the unrestricted fit alone and needs no refit,
which is why the summary reports it by default. It depends on how the
parameter is written: a nonlinear reparametrization of the coefficient
changes it, whereas the likelihood ratio, a difference of objective
values at the same two points, is not changed. The Wald statistic also
loses power against a distant alternative (the Hauck-Donner effect), the
curvature being read at \\\hat\beta\\ rather than under the null.

The gradient statistic needs the restricted fit and no matrix inversion:
it is the score contracted with the distance between the two fits.

## See also

[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
which the three restricted statistics read.
