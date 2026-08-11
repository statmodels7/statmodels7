# Estimate the Hyperparameters, by Whichever Route Each One Admits

Routes the twice differentiable hyperparameters to
[`outer_fit`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
and the rest to
[`statmod_path`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md).

## Usage

``` r
statmod_select(
  spec,
  design,
  blocks,
  hyper,
  inner_method,
  method,
  optimizer,
  beta,
  approx,
  maxit,
  tol,
  vb,
  data,
  weights,
  offsets
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- blocks:

  The blocks.

- hyper:

  The hyperparameters.

- inner_method:

  The inner method.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- optimizer:

  The optimizer for the differentiable hyperparameters.

- beta:

  The starting coefficients.

- approx:

  The approximation for the expected information.

- maxit, tol:

  The inner budget.

- vb:

  The verbosity.

- data, weights, offsets:

  What the fit was called on, which cross-validation needs in order to
  refit on part of it.

## Value

The list
[`outer_fit`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
returns.

## Details

The split is the same one that decides how the coefficients are fitted,
so a term whose penalty has a kink has both its coefficients and its
hyperparameters handled by methods that do not ask for a curvature it
does not have.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
