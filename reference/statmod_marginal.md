# The Marginal Criterion at Given Coefficients and Hyperparameters

The Laplace approximation to the log marginal likelihood, evaluated at
the penalized mode.

## Usage

``` r
statmod_marginal(
  spec,
  design,
  coef,
  hyper,
  method,
  approx = "bartlett",
  basis = NULL,
  ctx = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- approx:

  The approximation for the expected information.

- basis:

  The integrated subspace, from
  [`integrated_basis`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md).

## Value

A list with `value`, `loglik`, `penalty`, `logdet` and `q`, or `NULL`
where the determinant does not exist.
