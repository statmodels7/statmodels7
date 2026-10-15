# A Prediction-Error Criterion at Given Coefficients and Hyperparameters

Evaluates \\-2\ell + \kappa\tau\\ at one set of coefficients and
hyperparameters. The coefficients are taken as given and no fitting
happens here, so the criterion at the penalized mode is obtained by
passing the mode.

## Usage

``` r
statmod_pe(spec, design, coef, hyper, method, approx = "opg", active = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, refreshed at `coef` if any term needs it.

- coef:

  A named list of coefficient vectors, one per distribution parameter.

- hyper:

  The hyperparameters, per penalized term.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  of kind `"aic"` or `"bic"`, read for `hessian` and for \\\kappa\\.

- approx:

  How the expected information is approximated for a family with no
  closed form. Read only when `method@hessian` is `"expected"`.

- active:

  Which coefficients are away from a kink, as
  [`statmod_active()`](https://statmodels7.github.io/statmodels7/reference/statmod_active.md)
  reports them, or `NULL` where no penalty has one.

## Value

A list with the components below, or `NULL` when the penalized
information has no Cholesky factor or the trace is not finite:

- `value`:

  the criterion, \\-2\ell + \kappa\tau\\.

- `loglik`:

  the weighted log-likelihood at `coef`.

- `penalty`:

  the penalties' total value there. Reported for the caller; it does not
  enter `value`.

- `edf`:

  the trace \\\tau\\.

## Details

The pieces are assembled once each: the weighted log-likelihood, the
penalized information \\H + S\\ and the information \\H\\ alone, then
[`outer_tau()`](https://statmodels7.github.io/statmodels7/reference/outer_tau.md)
for the trace and
[`outer_k()`](https://statmodels7.github.io/statmodels7/reference/outer_k.md)
for \\\kappa\\.

The penalty's own value is reported alongside but does not enter the
criterion. A prediction-error criterion prices the fit's complexity
through \\\tau\\, not through the size of the penalty.

Where the penalty's Hessian is exactly zero on the active coordinates –
a lasso and nothing else penalizing them – the trace is
\\\mathrm{tr}(H\_{AA}^{-1}H\_{AA}) = \|A\|\\, the count of Zou, Hastie
and Tibshirani (2007), and it is returned as that count without forming
\\H\_{AA}\\. The condition is on the whole active block, so a ridge, a
random effect or a smooth among the active coordinates puts its
\\\lambda\\ into \\S\_{AA}\\ and the full trace is computed as before.

Where the active columns are linearly dependent \\H\_{AA}\\ is singular
and the trace does not exist; the degrees of freedom are then
\\\mathrm{rank}(X_A)\\ (Tibshirani and Taylor, 2012), which the count
overstates by the deficiency. It happens where a development carries its
own intercept beside a lasso over indicators that sum to it, as in
`nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))`. The count is therefore
taken only where
[`design_count_exact()`](https://statmodels7.github.io/statmodels7/reference/design_count_exact.md)
certifies that every active subset has full rank, and
[`active_rank()`](https://statmodels7.github.io/statmodels7/reference/active_rank.md)
answers everywhere else.

The certificate costs one pivoted decomposition per equation and is read
once where the design stands still, so a lasso path pays it at its first
evaluation only.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the criterion,
[`outer_tau()`](https://statmodels7.github.io/statmodels7/reference/outer_tau.md)
for the trace,
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
for the marginal criteria's counterpart.
