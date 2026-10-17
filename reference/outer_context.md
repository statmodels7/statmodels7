# One Evaluation Point, Shared

Creates the environment that holds the quantities that the marginal
criterion, its gradient and its Hessian all read at one \\(\beta,
\theta)\\: the linear predictors, the information, the penalty's
Hessian, their sum and its factorization.

## Usage

``` r
outer_context(spec, design, coef, hyper, approx = "opg")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- approx:

  The approximation for the expected information.

## Value

An environment carrying the point and, as they are requested, the
quantities derived from it.

## Details

The criterion, its gradient and its Hessian read the same information,
the same penalty Hessian and the same factorization at one point. The
context holds them, so that each is computed once.

The context is an environment, so an accessor fills it in place and
later readers find the quantity already there. Each quantity is computed
on first demand and never speculatively: the criterion alone does not
need an inverse, and a search that runs without an exact gradient does
not compute one.

Passing `NULL` wherever a context is accepted computes every quantity
afresh at the point given.

## See also

[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
