# Select the Hyperparameters of a Kinked Penalty Along a Path

Sweeps each of them over a grid of kink sizes, holding the others, and
keeps the setting the criterion prefers.

## Usage

``` r
statmod_path(
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
  offsets,
  rows,
  sweeps = 2L
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

- rows:

  Which hyperparameters to select.

- sweeps:

  How many cyclic passes.

## Value

The same list
[`outer_fit`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
returns.

## Details

The grid runs from the kink that empties the block down to `min_ratio`
of it, so the sweep goes from the sparsest fit towards the fullest and
every fit starts from the previous one's coefficients. Where the top of
the grid does not empty the block it is doubled until it does, the
starting value being computed at the coefficients in hand rather than at
a refitted null.

With several such hyperparameters the sweeps are cyclic, one coordinate
at a time, which is what keeps the cost linear in their number where a
full grid would be exponential in it.

Where the model also carries hyperparameters that are twice
differentiable, those are estimated by
[`outer_fit`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
inside each point of the path, so the two kinds are not mixed into one
search.

## See also

[`cv`](https://statmodels7.github.io/statmodels7/reference/cv.md),
[`path_rows`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
