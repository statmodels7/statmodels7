# Select the Hyperparameters of a Kinked Penalty Along a Path

Sweeps the hyperparameters of the kinked penalties over grids of kink
sizes and keeps the setting at which the criterion is best, by its
minimum or by the one-standard-error rule.

## Usage

``` r
statmod_path(
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
  rows,
  sweeps = 2L,
  nested_method = NULL
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
  that scores the points of the path.

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

- rows:

  The axes to select, as
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
  returns them.

- sweeps:

  The largest number of cyclic passes, `2` by default.

- nested_method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  that estimates the differentiable hyperparameters inside each point of
  the path, or `NULL`, in which case `method` is used.

## Value

A list with `par`, `hyper`, `value`, `criterion`, `converged`, `obj`,
`hist_blocks`, `hist_inner`, `hist_outer`, `iterations`, `evaluations`,
`exact_gradient` and `exact_hessian`: the fields of
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)'s
result that a path produces. `criterion` is the value of `method` at the
setting chosen, and `hist_outer` has one row for each point visited.

## Details

The grid runs from the kink that empties the block down to `min_ratio`
of it, so the sweep goes from the sparsest fit towards the fullest and
every fit starts from the previous one's coefficients. The top of the
grid is computed at the coefficients in hand and not at a refitted null;
where the fit at the top does not empty the block, the top is doubled
until it does.

A term that carries several of these hyperparameters has every
combination visited where it sets `search = "grid"` and one coordinate
at a time where it sets `"cyclic"`. Between terms the alternation is
cyclic either way, so that the cost is the product within a term and the
sum across terms. Each axis is built at the settings of the axes outside
it, which makes the elastic net's grid a family of \\\lambda\\ axes and
not a single one, and the axis swept by kink size is put innermost so
that the warm starts walk along it. The runs of a product grid are
independent and are distributed over worker processes where the
specification requests them.

A pass that would visit the points scored by the previous one is not
run. The top of the path is computed again at every pass, because the
rest of the model moves, and where it has not moved the grid is the one
already in hand. At most `sweeps` passes are made, and only where more
than one hyperparameter is selected.

Where the model also carries hyperparameters that are twice
differentiable, those are estimated by
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
inside each point of the path, by `nested_method`, so that the two kinds
are not mixed into one search. The unpenalized coefficients of a
dispersion are then held at the joint mode along the path when
`nested_method` is a marginal criterion, and are estimated on that
criterion once, at the setting chosen.

An error is signaled where `method` is
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) and
no point of a path is scored. A warning is signaled where the rule is
`"min"`, the axis is swept by kink size, and the minimum lies at the
dense end of the path with the criterion still falling there.

## See also

[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
[`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
