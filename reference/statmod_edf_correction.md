# The Smoothing-Parameter Correction to the Effective Degrees of Freedom

Returns the amount by which \\\mathrm{tr}\[(H+S)^{-1}H\]\\ understates
the complexity of a fit whose hyperparameters were themselves estimated.

## Usage

``` r
statmod_edf_correction(
  spec,
  coef,
  hyper,
  design,
  method,
  expected = TRUE,
  emethod = NULL,
  pinned = NULL
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

  Whether the information \\H\\ is the expected one. The argument is not
  read where the model carries a filter, whose information is observed.
  The caller passes the information that the fit used, so that
  \\V\_\beta\\, \\H\\ and the criterion's Hessian read one information.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

- pinned:

  The coordinates that the fit left out of the criterion's determinant,
  as it records them in `methods$pinned`.

## Value

A list with `total`, the scalar correction, `per`, one entry for each
parameter and term, and `n_hyper`, the number of hyperparameters that
were estimated. `per` uses only the diagonal of \\V\_\theta\\, so its
sum differs from `total` when the estimated hyperparameters are
correlated. The total and `n_hyper` are zero and `per` is empty where no
hyperparameter was estimated. A zero `total` beside a positive `n_hyper`
means that the curvature could not be read.

## Details

The ordinary effective degrees of freedom treat the smoothing parameters
as known, although they were chosen from the same data. Propagating
their uncertainty into the coefficients gives the corrected Bayesian
covariance

\$\$V' = V\_\beta + J V\_\theta J^\top, \qquad J =
\partial\hat\beta/\partial\theta,\$\$

and the corrected count is \\\mathrm{tr}(V' H)\\. \\J\\ comes from the
implicit function theorem at the penalized mode: differentiating
\\\partial(-\ell + \rho)/\partial\beta = 0\\ in the hyperparameter gives
\\(H+S)J_k = -\partial^2\rho/\partial\beta\\\partial\theta_k\\, whose
right-hand side is penalties7's `penalty_cross()`. \\V\_\theta\\ is the
inverse of the negative of the outer criterion's Hessian, which
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
returns; the negative is taken because the criterion is maximized. The
Hessian is that of the criterion that the fit maximized, so for
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
the determinant spans the range space of the penalty
([`integrated_basis()`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md)).

Everything is on the hyperparameter's link scale, which is the scale on
which the outer criterion is optimized and therefore the only scale on
which its Hessian is a variance.

The correction corresponds to the term of the effective degrees of
freedom `edf2` of mgcv that arises from the uncertainty of the smoothing
parameters (Wood, Pya and Safken, 2016). The scale of a Gaussian
response is a distribution parameter with its own coefficients here, so
its uncertainty is part of \\H\\ and needs no separate term.

The function returns zero, and `n_hyper` is zero, for every outer
criterion other than
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
including
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
used for a smooth penalty. A kinked penalty has no hyperparameter that
the outer criterion estimates, because
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
skips it, so there is nothing to propagate and the correction is zero.
The map from the hyperparameter to the penalized mode turns a corner
whenever a coefficient joins or leaves the active set, and a delta
method needs a derivative that does not exist there.

Where the criterion also estimates coefficients
([`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)),
the function returns
[`marginal_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/marginal_edf_correction.md).

The mode of a model that carries a filter moves in the joint vector of
the coefficients and the own parameters of the structural term, which is
the vector that the determinant of the criterion spans. A penalty over
those parameters is not a column of any design, so
[`hyper_mode_cross()`](https://statmodels7.github.io/statmodels7/reference/hyper_mode_cross.md)
omits it when it is called for the coefficients alone. The correction is
therefore computed on the joint vector.
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
takes the effective degrees of freedom of a structural model from
[`joint_smoother_diag()`](https://statmodels7.github.io/statmodels7/reference/joint_smoother_diag.md),
which reads the same two matrices, so the base count and the correction
are both computed on the joint vector. The two matrices are those that
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
and
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
build, so that this function and the criterion use one assembly. The
joint route reads the observed information, and `expected` and `emethod`
do not reach it. A model with no structural term of the filter shape
takes the coefficient route, because
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
returns `NULL` there.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
