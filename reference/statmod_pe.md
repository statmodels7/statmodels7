# A Prediction-Error Criterion at Given Coefficients and Hyperparameters

\\-2\ell + \kappa\tau\\ at the penalized mode.

## Usage

``` r
statmod_pe(spec, design, coef, hyper, method, approx = "bartlett")
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- approx:

  The approximation for the expected information.

## Value

A list with `value`, `loglik`, `penalty` and `edf`, or `NULL` where the
information is not invertible.

## See also

[`aic`](https://statmodels7.github.io/statmodels7/reference/aic.md)
