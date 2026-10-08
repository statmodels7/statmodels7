# The Derivative Row of an Equation Carrying a Filter

The derivative of the predictor of the equation carrying a score-driven
term, one row per observation, in every coordinate of the variance
matrix it moves with: the coefficients of every equation and the term's
free parameters on their unconstrained scale.

## Usage

``` r
structural_se_columns(spec, design, ep, p, coef)
```

## Arguments

- spec:

  The specification.

- design:

  Its design.

- ep:

  The predictors, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- p:

  The distribution parameter whose equation is being read.

- coef:

  The coefficients the predictors were evaluated at.

## Value

A list with `X` (the derivative rows) and `key` (the names of their
columns in the variance matrix), or `NULL` where the equation carries no
filter.

## Details

A filter's level is a recursion, not a column, so it has no row of a
design. Its derivative is the forward Jacobian of the recursion, which
[`filter_joint_jacobian()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_jacobian.md)
returns. Every equation's coefficients enter it, not only those of the
filter's own equation: the scores that drive the recursion are read at
the predictors of every equation, so a coefficient of the scale moves
the level of a filter in the mean. Leaving those columns out, on the
gaussian score-driven model of the Nile flow, gave a standard error of
the filtered mean up to 12 per cent too large.

A parameter that an intercept in the same equation holds is not
estimated and is not in that matrix, so it is not here either.

## See also

[`predict_se()`](https://statmodels7.github.io/statmodels7/reference/predict_se.md),
[`filter_joint_jacobian()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_jacobian.md)
