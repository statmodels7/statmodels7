# Estimate the Hyperparameters, by Whichever Route Each One Admits

Routes the twice differentiable hyperparameters to
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
and the rest to
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md).

## Usage

``` r
statmod_select(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  method,
  optimizer,
  beta,
  approx,
  maxit,
  tol,
  vb,
  data,
  sparse_method = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- blocks:

  The blocks.

- hyper:

  The hyperparameters.

- inner_optimizer:

  The inner method.

- method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  that estimates the twice differentiable hyperparameters, or `NULL`.

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

- data:

  The data the fit was called on, which cross-validation needs in order
  to refit on part of it; the weights and offsets are those held by
  `spec`.

- sparse_method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  that scores the path, or `NULL`, in which case `method` is used.

## Value

The list
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
or
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
returns, with the element `criterion_by`, the
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
whose value is `criterion`; the list
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
returns where no criterion runs.

## Details

The split is the same one that decides how the coefficients are fitted,
so a term whose penalty has a kink has both its coefficients and its
hyperparameters handled by methods that do not require a curvature that
the penalty does not have.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
