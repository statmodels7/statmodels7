# Fit a Separable Block by Coordinate Descent

Estimates one penalized block by cycling over its coefficients on the
working quadratic of its own equation while the other blocks are held
fixed. Where the equation has an unpenalized intercept that is not held,
the block is solved with the intercept profiled out, and the intercept
is then set to the value that goes with the block's coefficients.

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
  emethod,
  maxit = 100L,
  tol = 1e-08,
  prev_kink = NULL
)
```

## Arguments

- obj:

  The full objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it. The function reads the value at the point reached, not the
  gradient.

- beta:

  The current coefficients of all distribution parameters, stacked into
  one numeric vector (the vector that `obj$split()` divides by
  parameter). Every coefficient outside this block is held at its value,
  except the equation's unpenalized intercept when the block is
  centered.

- block:

  One entry of `statmod_blocks()$sparse`: the equation, the term, its
  column positions and its penalty.

- hyper:

  The hyperparameters, per penalized term, held fixed here.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it. Where a term recomputes its block with the coefficients,
  the block is read again at the current coefficients.

- expected:

  `TRUE` for the expected information in the working weights, `FALSE`
  for the observed one. Where the observed curvature of the equation is
  not positive at some observation, the expected information is used in
  its place.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

- maxit:

  The budget in weighted least squares iterations, each of which
  rebuilds the weights and runs the compiled sweeps to convergence.

- tol:

  The stopping tolerance on the largest absolute change in the block's
  coefficients over a weighted least squares iteration. The compiled
  sweeps use the same tolerance on the largest absolute change of a
  sweep.

- prev_kink:

  The size of the kink at the previous point of a path, a single number,
  or `NULL` to cycle over every coordinate. Only the strong rule reads
  it.

## Value

A list shaped like
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)'s:
`par` (the whole stacked coefficient vector, updated in this block and,
where the block is centered, in the equation's intercept), `value` (the
objective), `converged`, `iterations` (the total number of kernel
sweeps) and `method` (`"coordinate descent"`). `converged` is `TRUE`
when the iterations settled or when the objective did not fall for any
length of step. The value is `NULL` where the route does not apply: the
penalty describes no proximal table at the steps that this block's
curvature produces, or the working weights are unusable.

## Comparison with the proximal method

Each iteration of a proximal gradient method evaluates the objective,
its gradient and the penalty's operator over every parameter of the
distribution. A coordinate descent reads only the block's own columns
and the running residual.

## The working quadratic

With \\\eta\\ the equation's linear predictor, \\s_i\\ the score in it
and \\h_i\\ the information, \\-\ell\\ is \\\frac12\sum_i h_i(z_i -
\eta_i)^2\\ up to a constant with \\z = \eta + s/h\\, which is the
weighted least squares problem of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
restricted to one equation. The other columns of that equation enter as
an offset. For a Gaussian response with an identity link the quadratic
is exact and the first weighted least squares iteration gives the
answer; otherwise the weights are rebuilt and the sweeps repeated.

## Centering

Where the equation has an unpenalized intercept that is not held (see
[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)),
the compiled sweeps solve the block with the intercept profiled out: the
columns and the working response enter at their weighted means,
implicitly, so a sparse block stays sparse. After the sweeps the
intercept is set to the weighted mean of the working response net of the
block, so the step is a joint step in the block and the intercept. A
held intercept, or an equation with no intercept, leaves the block
uncentered.

## The penalty arrives as a table

The coordinate update is the penalty's own proximal operator at the step
\\1/v_j\\, with \\v_j = \sum_i w_i x\_{ij}^2\\ (the centered sum of
squares where the block is centered), and \\v_j\\ does not move while
the working weights are held. The whole table is therefore built once
per weighted least squares iteration by
[`penalties7::penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
and the compiled sweeps read it, so the kernel does not refer to any
family by name, and a penalty that describes its operator uses the
compiled route without an edit to the kernel.

## Screening

In passing from one point of a path to the next, a coordinate can be
discarded when the gradient that it had at the previous point is below
\\2s_k - s\_{k-1}\\, with \\s\\ the size of the kink: the sequential
strong rule of Tibshirani and others (2012), which assumes the gradient
moves at most as fast as the threshold does. That assumption is not a
theorem, so the rule can discard a coordinate that belongs in the fit.
The answer is exact because of the check afterwards: the gradient is
read over every column at the point reached, any discarded coordinate
whose gradient exceeds the kink is put back, and the fit is repeated.
Without the check the answer would not be exact.

## Which update

The gradient is kept either as a residual, at \\O(n)\\ a visit, or as
itself through \\g_j = (X'Wz)\_j - \sum_k (X'WX)\_{jk}\beta_k\\, at
\\O(m)\\ a change with the Gram columns cached as coordinates come
alive. The second costs less when \\n\\ is large next to the number of
live coordinates and uses more memory, so
[`coord_covariance()`](https://statmodels7.github.io/statmodels7/reference/coord_covariance.md)
computes the choice from the two sizes.

## Acceptance

A step is accepted where the penalized objective falls, and is halved
towards the current point otherwise. A SCAD or MCP penalty that carries
a curvature (see
[`coord_table_penalty()`](https://statmodels7.github.io/statmodels7/reference/coord_table_penalty.md))
takes every step whole.

## References

Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. *Journal of
Statistical Software* 33(1), 1–22.

Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor,
J. and Tibshirani, R. J. (2012). Strong rules for discarding predictors
in lasso-type problems. *Journal of the Royal Statistical Society,
Series B* 74(2), 245–266.

## See also

[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md),
[`penalties7::penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
