# Directions Along Which the Estimate Does Not Exist

Searches the coefficients of a fit for a direction of recession, along
which the penalized log-likelihood keeps increasing, or stays level, as
the linear predictors move without bound. Such a direction exists when
the maximum likelihood estimate does not exist (separation in a binary
or binomial response, a level of a factor with all its counts at zero, a
group without zeros in the equation of a zero-inflation probability),
and the fit then stops at a finite point that depends only on its
stopping rule.

## Usage

``` r
separation_diagnosis(
  fit,
  steps = 2^(0:6),
  tol = 1e-08,
  flat = 1e-06,
  max_eigen = 4L
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- steps:

  The multiples of the scaled direction at which the objective is
  evaluated, in units of the linear predictor.

- tol:

  The largest rise of the objective, relative to \\\max(1,
  \|F(\hat\beta)\|)\\, that is still read as level.

- flat:

  The largest ratio of an eigenvalue to the largest one at which an
  eigenvector of \\S\\ against \\M\\ is probed beside the Newton step.

- max_eigen:

  The largest number of eigenvectors that are probed.

## Value

A list with `found` (logical, a direction of recession), `descent`
(logical, a direction of descent), `directions` (one entry per direction
found, each with `kind`, `"recession"` or `"descent"`; `coefficients`,
the named components of the direction over the coefficients that it
moves, scaled so that the largest change in a linear predictor is one;
`free`, the names of those coefficients that carry no smooth penalty;
`equations`; `n_moved`, the number of observations whose linear
predictor moves by more than 1e-3 in each equation; `change`, the change
of the objective at each step; `source`, `"newton"` or `"eigen"`; and
`curvature`, the curvature along the direction against \\M\\ as a ratio
to the largest eigenvalue), `newton_eta`, the largest change that the
Newton step makes in a linear predictor, `min_eigen`, the smallest ratio
of an eigenvalue of \\S\\ to the largest one, `n_eval`, the number of
evaluations of the objective, and `reason`, a character string where the
check could not be made.

## Details

The objective is \\F\\ of
[`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md),
at the fitted hyperparameters, with gradient \\g\\ and Hessian \\K\\.
The coefficients are split into those that carry no smooth penalty
(\\u\\) and those that do (\\v\\), and the latter are profiled out
through the Schur complement \\S = K\_{uu} - K\_{uv} K\_{vv}^{-1}
K\_{vu}\\, so that the eigenvalue problem has the size of the
unpenalized coefficients whatever the number of random effects. A
penalized coefficient moves along each candidate as the profile
requires, \\v = -K\_{vv}^{-1} K\_{vu} u\\, and only the null space of a
penalty can carry a direction of recession.

The candidates are the Newton step at the fitted point, \\\delta =
-K^{-1} g\\, with each eigenvalue of \\S\\ floored at \\10^{-14}\\ times
the largest; where coefficients are estimated on the outer criterion, a
second Newton step that keeps their gradient in \\F\\ (see below); and
both signs of the eigenvectors of \\S\\ against the metric of the linear
predictors, \\M = X_u^\top X_u\\ over the equations, whose eigenvalues
fall below `flat` times the largest one. At a point where the fit has
drifted along a direction of recession the gradient and the curvature
along it vanish together, so the Newton step remains of order one in the
linear predictor while the decrease that it predicts is below the
stopping rule.

Each candidate \\d\\ is scaled so that the largest change it makes in a
linear predictor is one, and \\F\\ is evaluated at \\\hat\beta + t d\\
for every `t` in `steps`. A candidate is a direction of recession when
\\F\\ does not rise by more than `tol` times \\\max(1,
\|F(\hat\beta)\|)\\ at any of those points, nor from one point to the
next, and the whole decrease stays within
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).
It is a direction of descent when \\F\\ falls by more than
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
at one of them, which means that the point is not a minimum at all (a
scale whose group has a constant response, where the log-likelihood is
unbounded, is an example). At an ordinary optimum a change of one in a
linear predictor already raises \\F\\.

Coefficients held at the kink of a penalty, those of a block that is a
working linearization
([`frozen_block()`](https://statmodels7.github.io/statmodels7/reference/frozen_block.md))
and aliased coefficients are held fixed. Coefficients estimated on the
outer criterion move with the others, since a direction of recession of
\\F\\ is one of that criterion too when no random effect integrates over
them, but their gradient in \\F\\ is set to zero, because they sit at
the root of the outer criterion and not of \\F\\. For the same reason a
direction of descent is read on the change of \\F\\ after removing the
linear term of those coordinates, \\t\\g_m^\top d\\, where \\g_m\\ is
their gradient in \\F\\: under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
the scale of a gaussian is not at the minimum of \\F\\, and without the
correction the difference between the two criteria reads as a direction
of descent. A direction of recession is read on \\F\\ itself.
