# Fit a Separable Block by Coordinate Descent

Estimates one penalized block by cycling over its coefficients on the
working quadratic of its own equation, the other blocks held fixed.

## Usage

``` r
coord_fit(
  obj,
  beta,
  block,
  hyper,
  spec,
  design,
  expected,
  approx,
  maxit = 100L,
  tol = 1e-08
)
```

## Arguments

- obj:

  The full objective.

- beta:

  The current stacked coefficients.

- block:

  One entry of `statmod_blocks()$sparse`.

- hyper:

  The hyperparameters.

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- expected:

  Whether the information is the expected one.

- approx:

  How the expected information is approximated.

- maxit:

  The iteration budget.

- tol:

  The stopping tolerance.

## Value

A list shaped like
[`sparse_fit`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)'s,
or `NULL` where the route does not apply.

## Details

**Why not the proximal method.** A proximal gradient step reads the
whole model: measured on 200 observations and 20 columns, one block fit
made 88 evaluations of the objective, 75 of the gradient and 83 of the
operator, each over every parameter of the distribution, and closed in
36 iterations at 0.17 seconds. A coordinate descent reads the block's
own columns and the running residual instead and closes in six sweeps.

**The working quadratic.** With \\\eta\\ the equation's linear
predictor, \\s_i\\ the score in it and \\h_i\\ the information,
\\-\ell\\ is \\\frac12\sum_i h_i(z_i - \eta_i)^2\\ up to a constant with
\\z = \eta + s/h\\, which is the weighted least squares problem of
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
restricted to one equation. The other columns of that equation enter as
an offset. For a Gaussian response with an identity link the quadratic
is exact and one pass is the answer; otherwise the weights are rebuilt
and the sweeps repeated.

**The penalty arrives as a table.** The coordinate update is the
penalty's own proximal operator at the step \\1/v_j\\, with \\v_j =
\sum_i w_i x\_{ij}^2\\, and \\v_j\\ does not move while the working
weights are held. The whole table is therefore built once per weighted
least squares iteration by
[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
and the compiled sweeps read it, so the kernel names no family and a
penalty that describes its operator gets the compiled route without an
edit here.

## See also

[`sparse_fit`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md),
[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
