# What the Derivative of a Continued Filter Starts From

The derivatives of the level and of the score at the observed rows, and
of the term's parameters, in the coordinates the variance matrix of the
fit is written in: the coefficients of every equation followed by the
term's free parameters on their unconstrained scale.

## Usage

``` r
continued_deriv_inputs(ospec, odesign, coef, f, ost, tm)
```

## Arguments

- ospec, odesign:

  The fit's specification and design.

- coef:

  The coefficients.

- f:

  The filter object at the observed rows.

- ost:

  The structural state.

- tm:

  The term.

## Value

A list with `deriv` (as
[`modelterms7::term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.html)
takes it), `col` (the columns of the filter equation's coefficients) and
`key` (the names of the columns in the variance matrix), or `NULL` where
the fit carries no filter.

## Details

The coordinates are the ones
[`structural_se_columns()`](https://statmodels7.github.io/statmodels7/reference/structural_se_columns.md)
uses for a prediction at the observed rows, so a forecast and a fitted
value are read against the same rows and columns of the variance. The
derivative of the filtered predictor at an observed row is the forward
Jacobian of
[`filter_joint_jacobian()`](https://statmodels7.github.io/statmodels7/reference/filter_joint_jacobian.md);
the level's derivative is that minus the equation's design row. The
score drives the recursion and depends on the predictors of every
equation, so its derivative is \\\sum_b \ell\_{pb,t} V\_{b,t}\\, with
\\V\_{b,t}\\ the derivative row of equation \\b\\ and \\\ell\_{pb,t}\\
the family's second derivative on the link scale. A parameter's
derivative in its own unconstrained coordinate is its link's, and one
for a coefficient of a development.

## See also

[`statmod_eta_continued()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta_continued.md),
[`structural_se_columns()`](https://statmodels7.github.io/statmodels7/reference/structural_se_columns.md)
