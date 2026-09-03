# A Prediction-Error Criterion at Given Coefficients and Hyperparameters

Evaluates \\-2\ell + \kappa\tau\\ at one set of coefficients and
hyperparameters. The coefficients are taken as given and no fitting
happens here, so a caller who wants the criterion at the penalized mode
has to pass the mode.

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

A list of four, or `NULL` when the penalized information has no Cholesky
factor:

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

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the criterion,
[`outer_tau()`](https://statmodels7.github.io/statmodels7/reference/outer_tau.md)
for the trace,
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
for the marginal criteria's counterpart.
