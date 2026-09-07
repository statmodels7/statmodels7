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
  mode_error = FALSE
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
  ratio and the gradient statistic need none, and the two loops that
  call this in quantity –
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
  one row at a time, and
  [`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
  five to seven times per interval – do not read it.

## Value

A list with `test`, `statistic`, `df`, `p.value`, `converged` – `NA` for
the Wald statistic and the restricted fit's flag for the other three –
and `mode_error`, `NA` unless it was asked for.

## Details

All four read the PENALIZED objective, which is the function the fit
maximizes. Writing \\\ell_p(\beta) = \ell(\beta) - \rho(\beta;\theta)\\
for it, \\\hat\beta\\ for the unrestricted estimates, \\\tilde\beta\\
for those with \\\beta_j\\ held at \\b\\, \\U = \partial\ell_p\\ for its
score and \\K\\ for the penalized information, the four are \$\$W =
(\hat\beta_j - b)^2 / \widehat{\mathrm{Var}}(\hat\beta_j), \qquad LR =
2\\\ell_p(\hat\beta) - \ell_p(\tilde\beta)\\,\$\$ \$\$S =
U_j(\tilde\beta)^2 \\ \[K(\tilde\beta)^{-1}\]\_{jj}, \qquad T =
U_j(\tilde\beta)\\(\hat\beta_j - b),\$\$ each compared with a
\\\chi^2_1\\. Where the model carries no penalty \\\ell_p\\ IS the
log-likelihood and these are the classical four.

The score statistic's general form is \\U'K^{-1}U\\, which reduces to
the product above because every other component of \\U\\ vanishes at the
restricted maximum – and that is worth knowing rather than assuming,
since it is what says the restricted fit converged.

## Why the penalized objective and not the likelihood

Reading the likelihood alone would make the four incoherent with one
another and the second of them not a ratio at all. \\\tilde\beta\\
maximizes \\\ell_p\\ under the restriction and not \\\ell\\, so
\\\ell(\tilde\beta)\\ can exceed \\\ell(\hat\beta)\\ and the difference
come out NEGATIVE – measured on a smooth, holding `s(z).z3` at two
values either side of its estimate gives -0.0379 and -0.1075 where
\\\ell_p\\ gives +0.0054 and +0.1355. The other three were already
reading the penalized objective, \\U\\ and \\K\\ being its gradient and
its curvature and the Wald variance being \\(H+S)^{-1}\\, so one table
would have carried two different tests.

What it means for a SHRUNK coordinate is the reading
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
already gives it under `type = "bayesian"`: the penalty is a log-prior,
so \\\ell_p\\ is a log posterior at the hyperparameters in force and the
interval by inversion is the set the posterior does not reject,
conditional on them. Nothing new is claimed; the claim the Wald interval
on such a row already carries is extended to the other three. Measured,
the two agree as they must: on a smooth at \\n = 300\\ the inverted
interval and the Wald one differ by 8.7e-04 on `s(z).z1`, 1.1e-03 on an
unpenalized covariate and 2.5e-03 on the smooth's linear column, which
is second-order agreement, and what the inversion adds is the asymmetry.

## What separates them

Wald reads the unrestricted fit alone and needs no refit, which is why
it is the one every summary prints and also the one that is NOT
invariant to how the parameter is written: a reparametrization moves it
and moves none of the other three. It is also the one that loses power
against a distant alternative, the Hauck-Donner effect, the curvature
being read at \\\hat\beta\\ rather than under the null.

The gradient statistic needs the restricted fit but NO matrix inversion
at all, which is its whole point: it is the score contracted with the
distance between the two fits.

## See also

[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
which the three restricted statistics read.
