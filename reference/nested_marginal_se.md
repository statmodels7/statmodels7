# The Standard Error of a Marginal Parameter With Nested Effects

The delta method of \\\bar\theta = \sum_k w_k h^{-1}(\eta(\beta,
u_k))\\, whose gradient in the coefficients is \\\sum_k w_k\\
(h^{-1})'(\eta_k)\\ J_k\\ with \\J_k\\ the design at the coefficients of
node \\k\\.

## Usage

``` r
nested_marginal_se(fit, pr, coef_at, eta_at, w, p, fitv, level = 0.95, ...)
```

## Arguments

- fit:

  The fitted model.

- pr:

  What
  [`nested_prepare()`](https://statmodels7.github.io/statmodels7/reference/nested_prepare.md)
  returns.

- coef_at:

  A function of the node index returning the coefficients there.

- eta_at:

  A function of the node index returning the predictors there.

- w:

  The nodes' weights.

- p:

  The parameter.

- fitv:

  The marginal parameter.

- level:

  The confidence level.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A data frame with `fit`, `se`, `lower` and `upper`; the interval is the
delta method's on the link scale, carried through the link.
