# Which Coefficients Are Not Sitting at a Kink

A logical vector over the stacked coefficients, `FALSE` where one lies
at a point its penalty is not differentiable at.

## Usage

``` r
statmod_active(spec, blocks, beta, hyper, tol = 1e-08)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The blocks, as
  [`statmod_blocks`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns them.

- beta:

  The stacked coefficients.

- hyper:

  The hyperparameters.

- tol:

  How close to a kink counts as at it.

## Value

A logical vector as long as `beta`.

## Details

The kink locations come from
[`penalty_kinks`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
read at the hyperparameters in force, and a coefficient is inactive when
it sits at one of them. Everything outside a kinked block is active,
having no kink to sit at.

The map of a kinked penalty is the identity here: a separable penalty
under a general map is the generalized-lasso problem, which
[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.html)
rejects, so a block that reaches this function penalizes its
coefficients one at a time and a kink of the penalty is a kink in a
coefficient.

## See also

[`outer_tau`](https://statmodels7.github.io/statmodels7/reference/outer_tau.md)
