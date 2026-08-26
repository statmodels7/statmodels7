# The Columns a Structural Term Adds to a Derivative Row

The derivative of one equation's predictor in the free parameters of the
score-driven term sitting in it, one row per observation, on the scale
the variance matrix is indexed by.

## Usage

``` r
structural_se_columns(spec, design, ep, p, X)
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

- X:

  The equation's design rows.

## Value

A list with `X`, `J` and `key`, or `NULL` where the equation carries no
filter.

## Details

A filter's level is a recursion, not a column, so it has no row of a
design. What it has is the derivative the recursion propagates beside
the state, which
[`modelterms7::term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.html)
returns on the parameter scale; multiplying by each parameter's own
\\h'(\zeta_j)\\ carries it to the unconstrained scale the joint matrix
is written in, which is the chain the exact gradient already uses.

A parameter an intercept in the same equation holds is not estimated and
is not in that matrix, so it is not here either.

The design's own rows are corrected at the same time, through
[`modelterms7::term_static_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.html):
a coefficient of this equation moves the level as well as the static
part, because the scores driving the recursion are read at the predictor
the recursion produces. Measured on a score-driven mean with one
covariate beside it, leaving that out understates the standard error by
about a quarter.

## See also

[`predict_se()`](https://statmodels7.github.io/statmodels7/reference/predict_se.md),
[`statmod_filter_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
