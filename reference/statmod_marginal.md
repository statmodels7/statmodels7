# The Marginal Criterion at Given Coefficients and Hyperparameters

The Laplace approximation to the log marginal likelihood, evaluated at
the penalized mode: the log-likelihood minus the penalty, plus
\\(q/2)\log 2\pi\\, minus half the log determinant of the penalized
information over the integrated subspace. A coordinate with a non-finite
curvature is pinned to the identity (see
[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md)),
which lowers \\q\\ by one. The determinant is read through
[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md).

## Usage

``` r
statmod_marginal(
  spec,
  design,
  coef,
  hyper,
  method,
  approx = "opg",
  basis = NULL,
  ctx = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- approx:

  The approximation for the expected information.

- basis:

  The integrated subspace, from
  [`integrated_basis()`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md).

- ctx:

  An evaluation context from
  [`outer_context()`](https://statmodels7.github.io/statmodels7/reference/outer_context.md)
  built at the same coefficients and hyperparameters, or `NULL` to
  assemble the information and the penalty afresh.

## Value

A list with `value`, `loglik`, `penalty`, `logdet` and `q`, or `NULL`
where the determinant does not exist. Where the integrated subspace is
empty, `value` is the profile log-likelihood minus the penalty, with
`logdet` and `q` equal to zero.
