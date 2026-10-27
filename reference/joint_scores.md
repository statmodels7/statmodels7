# The Scores of the Observations in the Joint Vector of Unknowns

Returns the matrix whose row \\t\\ is the derivative of the
log-likelihood contribution of observation \\t\\, weighted by its prior
weight, in the coefficients of every equation followed by the own
parameters of a structural term on the unconstrained scale.

## Usage

``` r
joint_scores(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

## Value

A matrix with `spec@n_obs` rows and a column for each unknown of the
joint information.

## Details

Without a structural term the row is the family's score on the link
scale times the design row of each equation. With a filter it is the sum
over the equations of the score times the derivative of the predictor in
the unknowns, and the derivative of the filtered equation's predictor is
the forward Jacobian that
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
returns. With a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term it is the derivative of the increment of the forward recursion, see
[`regime_scores()`](https://statmodels7.github.io/statmodels7/reference/regime_scores.md).
The columns of a level that an intercept holds are dropped, as they are
from the information.
