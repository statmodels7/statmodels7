# The Smoothing-Parameter Correction to the Effective Degrees of Freedom

The amount by which \\\mathrm{tr}\[(H+S)^{-1}H\]\\ understates the
complexity of a fit whose hyperparameters were themselves estimated.

## Usage

``` r
statmod_edf_correction(
  spec,
  coef,
  hyper,
  design,
  method,
  expected = TRUE,
  approx = "opg"
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- design:

  The design.

- method:

  The outer method that estimated them, or `NULL`.

- expected:

  Whether the information is the expected one.

- approx:

  The approximation for the expected information.

## Value

A list with `total`, the scalar correction, `per`, one entry per penalty
key, and `n_hyper`, how many hyperparameters were estimated. Zero
throughout where none was; a zero `total` beside a positive `n_hyper`
means the curvature could not be read, which is what a shared
hyperparameter leaves, and a caller reporting to a reader has to tell
the two apart.

## Details

The ordinary effective degrees of freedom read the smoothing parameters
as though they were known, and they were not: they were chosen from the
same data. Propagating their uncertainty into the coefficients gives the
corrected Bayesian covariance

\$\$V' = V\_\beta + J V\_\theta J^\top, \qquad J =
\partial\hat\beta/\partial\theta,\$\$

and the corrected count is \\\mathrm{tr}(V' H)\\. \\J\\ comes from the
implicit function theorem at the penalized mode: differentiating
\\\partial(-\ell + \rho)/\partial\beta = 0\\ in the hyperparameter gives
\\(H+S)J_k = -\partial^2\rho/\partial\beta\\\partial\theta_k\\, whose
right-hand side is penalties7's `penalty_cross()`. \\V\_\theta\\ is the
inverse of the outer criterion's own Hessian, which
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
returns, negated because that criterion is a maximand.

Everything is on the hyperparameter's link scale, which is where the
outer criterion optimizes and therefore the only scale on which its
Hessian is a variance.

**Against mgcv.** This is the first of the two terms mgcv sums into
`edf2`, and it agrees with mgcv's to about 3e-4 on a univariate smooth
once the difference between the two bases is allowed for. mgcv adds a
second term for the Gaussian scale, which it profiles out of the fit and
whose uncertainty it must therefore add back; here every distribution
parameter carries its own equation and its own coefficients, so that
uncertainty is already inside \\H\\. The residual difference is measured
at 0.16, 0.10 and 0.05 effective parameters at n = 200, 400 and 2000,
falling with the sample size.

**Where it does not apply.** A kinked penalty has no hyperparameter the
outer criterion estimates,
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
skipping it, so there is nothing to propagate and the correction is
zero. That is not an approximation: the map from the hyperparameter to
the penalized mode turns a corner whenever a coefficient joins or leaves
the active set, and a delta method needs a derivative that does not
exist there.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
