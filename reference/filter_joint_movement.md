# How the Mode Moves Where an Unpenalized Filter Moves With It

The pieces
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
reads on the joint vector of coefficients and a filter's own parameters
where no penalty covers those parameters: a solve with the joint
penalized curvature, which is how the mode moves, and the derivative of
every equation's predictor over that vector, which is what the
determinant's movement is read along.

## Usage

``` r
filter_joint_movement(spec, design, coef, S, Dm, total, D3)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design as it arrived, carrying the structural state.

- coef:

  The coefficients at the penalized mode.

- S:

  The penalty's Hessian over the stacked coefficients.

- Dm:

  The mode's own correction for a block that moves, from
  [`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md).

- total:

  The number of stacked coefficients.

- D3:

  The family's third derivative on the link scale at the fitted
  predictors, which the recursion reads where the memo misses.

## Value

A list with `rows`, one matrix per distribution parameter over the
estimated coordinates, and `solve`, a function of a vector over the
coefficients returning the joint solve over every estimated coordinate;
or `NULL` where the joint curvature cannot be formed or solved.

## Details

The mode is where the penalized objective's gradient in the coefficients
AND in the filter's parameters vanishes, so differentiating that
condition in a hyperparameter gives \$\$(J + S)\\v =
-\partial^2\rho/\partial u\\\partial\theta,\$\$ with \\J\\ the joint
observed information of
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
and \\S\\ the penalty's Hessian placed on the coefficients, the only
coordinates it covers. A block that moves with its coefficients adds its
own second derivative on those coordinates, as it does on the
coefficient route.

The rows are the static design of every equation, with the equation
carrying the filter replaced by the forward Jacobian of the recursion,
read from
[`filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/filter_curvature.md)
at the same memo slot the joint information fills, so the recursion runs
once at a point. A held level is dropped from both, exactly as the
information drops it.

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/filter_curvature.md),
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
