# The Exact Derivatives of a Prediction-Error Criterion

The gradient and, when asked, the Hessian of \\-2\ell + \kappa\tau\\
over the free scale of the hyperparameters.

## Usage

``` r
statmod_pe_derivs(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  order = 1L,
  approx = "opg"
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The outer index.

- order:

  `1` for the gradient alone, `2` for both.

- approx:

  How the expected information is approximated for a family with no
  closed form. Read only when `method@hessian` is `"expected"`.

## Value

A list with `grad` and, at order 2, `hess`; or `NULL` where the
information is not invertible.

## Details

Writing \\P = J^{-1}\\, \\A_m = \partial J/\partial\theta_m = S_m +
T\[\hat\beta_m\]\\ and \\B_m = \partial H/\partial\theta_m =
T\[\hat\beta_m\]\\, \$\$\tau_m = \mathrm{tr}(PB_m) -
\mathrm{tr}(PA_mPH),\$\$ and the log-likelihood contributes
\\(\partial\rho/\partial\beta)' \hat\beta_m\\. Differentiating once more
brings in \\A\_{ml}\\, \\B\_{ml}\\ and \\\hat\beta\_{ml}\\, which are
the quantities
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
already assembles.

Where `method@hessian` is `"expected"`, \\\tau\\ reads the expected
information \\H_E\\, so its derivative contracts the derivative of
\\H_E\\ in the predictors, read from
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
under its own key, while the mode still moves by the observed penalized
information. That route has no Hessian and is refused at order 2.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
